!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use cooling_module, ONLY: twopi, kb, mh
  use incgamma
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0 || NVAR>NDIM+2+NENER
  integer::ivar
#endif
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  integer::i,ilun
  real(dp)::currad,xc,yc,zc,rho0g,rho0dm,rho_cloud,c_s2,PhiR,P_wind,P_cloud,nH
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_prs
  real(dp),save::mu_cloud,mu_wind,Phi0
  real(dp),save::r_max=0.0,velinit=0.0,rtinit=0.0
  real(dp)::tinit
  character(len=256)::fileloc
  logical::file_exists = .false.
  logical,save::firstcall=.true.
#ifdef EINASTO
  real(dp),save::gamma3n,ein_M
#endif

#ifdef SIMPLE_IC
  call region_condinit(x,q,dx,nn)
  q(1:nn,imetal) = Z_cloud*0.02
  q(1:nn,imetal+1) = 0.0	!delayed cooling
  q(1:nn,imetal+2) = 1.0	!tracer
#else
  ! User-defined initial conditions

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
  scale_prs = scale_d * scale_v**2

  xc = x1_c*boxlen
  yc = x2_c*boxlen
  zc = x3_c*boxlen
  if (firstcall) then
     nH = n0g*scale_nH
     call GetMuFromTemperature(T_cloud,nH,mu_cloud)
     nH = ndens_wind*scale_nH
     call GetMuFromTemperature(T_wind,nH,mu_wind)

     if (vel_wind > 0.0) then
        fileloc=trim(output_dir)//trim(orbitfile)
        inquire(file=fileloc,exist=file_exists)
        if(file_exists) then
           open(newunit=ilun, file=fileloc)
           read(ilun,*)tinit,velinit
           close(ilun)
        else
           write(*,*)"File ",trim(orbitfile)," missing!"
        endif
     else
        velinit = vel_wind
     endif

     if (evolve_rtidal)then
        fileloc=trim(output_dir)//trim(rtidalfile)
        inquire(file=fileloc,exist=file_exists)
        if(file_exists) then
           open(newunit=ilun, file=fileloc)
           read(ilun,*)tinit,rtinit
           close(ilun)
        else
           write(*,*)"File ",trim(rtidalfile)," missing!"
        endif
     endif

#ifdef EINASTO
     gamma3n = cmpgamma(3d0*ein_n)
     ein_M = 2d0*twopi*rhodm0*R_s**3*ein_n*gamma3n
     Phi0 = -ein_M*cmpgamma(2d0*ein_n)/(R_s*gamma3n)
#else
     ! NFW
     Phi0 = -2d0*twopi*rhodm0*R_s**2
#endif
     if(M_plummer > 0.0)Phi0 = Phi0 - M_plummer/r_plummer
     firstcall = .false.
  endif
  
  rho0g=n0g*mu_cloud
  c_s2 = kb*T_cloud/(mu_cloud*mh)/scale_v**2 !square of isothermal sound speed in cloud centre
  P_wind = ndens_wind*T_wind/scale_T2

  do i=1,nn
#if NDIM==3
    currad = dsqrt((x(i,1)-xc)**2 + (x(i,2)-yc)**2 + (x(i,3)-zc)**2)
#else
    currad = dsqrt((x(i,1)-xc)**2 + (x(i,2)-yc)**2)
#endif

#ifdef EINASTO
    PhiR = -ein_M*(R_s*gamma3n + currad*gammainc2n((currad/R_s)**(1d0/ein_n)) - R_s*gammainc3n((currad/R_s)**(1d0/ein_n)))/(R_s*currad*gamma3n)
#else
    ! NFW
    PhiR = Phi0*R_s*dlog(1+currad/R_s)/currad
#endif
    if (M_plummer > 0.0)PhiR = PhiR - M_plummer/sqrt(r_plummer**2 + currad**2)
    if (inner_dens > 0.0) then
       if (currad < r_inner) then
          rho_cloud = inner_dens*dexp(-inner_slope*currad)
       else
          rho_cloud = outer_dens*dexp(-outer_slope*currad)
       endif
    else
       rho_cloud = rho0g*dexp(-(PhiR-Phi0)/c_s2)
    endif
    P_cloud = (rho_cloud*T_cloud/mu_cloud)/scale_T2
    if ((evolve_rtidal .and. (currad < rtinit)) .or. ((.not.(evolve_rtidal)) .and. (P_cloud < P_wind)))then
      q(i,1) = ndens_wind*mu_wind
      q(i,ndim+2) = P_wind
      q(i,2) = 0.0      !x-velocity
      q(i,3) = velinit ! vel_wind*1.e5/scale_v !y-velocity (given in km/s)
#if NDIM==3
      q(i,4) = 0.0      !z-velocity
#endif
      q(i,imetal) = Z_wind*0.02
      q(i,imetal+1) = 0.0	!delayed cooling
      q(i,imetal+2) = 0.0	!tracer
!      q(i,imetal+3) = 0.0     !sf gas tracer
    else
      q(i,1) = rho_cloud
      q(i,ndim+2) = P_cloud
      q(i,2) = 0.0      !x-velocity
      q(i,3) = 0.0      !y-velocity
#if NDIM==3
      q(i,4) = 0.0      !z-velocity
#endif
      q(i,imetal) = Z_cloud*0.02 !metallicity converted from relative to solar to absolute assuming Z_sol=0.02 as hardcoded in other parts of RAMSES
      q(i,imetal+1) = 0.0	!delayed cooling
      q(i,imetal+2) = 1.0 !tracer
!      if (currad < Rad_cloud) then
!        q(i,imetal+3) = 1.0     !sf gas tracer
!      else
!        q(i,imetal+3) = 0.0     !sf gas tracer
!      endif
!      if (currad > r_max)r_max=currad
    endif
!    write(*,*) "i ", i, " r ", currad, " rho: ", q(i,1), " P:",q(i,ndim+2), "Phi:",PhiR,"gammainc2n:",gammainc2n((currad/R_s)**(1d0/ein_n)),"gammainc3n",gammainc3n((currad/R_s)**(1d0/ein_n))!" vy: ", q(i,3), " P: ", q(i,ndim+2), " T ", (q(i,ndim+2)*mu/q(i,1))*scale_t2, " x1: ", x1_c
  enddo
!write(*,*)'r_max:',r_max, 'mu_wind:',mu_wind,'mu_cloud:',mu_cloud
!endif
#endif
  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,ndim+2+ivar)=q(1:nn,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit

subroutine GetMuFromTemperature(T,nH,mu)
!Note: T is assumed here to be actual temperature in Kelvin, NOT T/mu
  use amr_parameters, ONLY: dp, aexp
  use cooling_module, ONLY: set_rates, cmp_chem_eq
  implicit none
  real(dp)::T,nH,mu
  real(dp)::mu_old,err_mu,mu_left,mu_right,n_TOT
  real(dp),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(dp),dimension(1:6) :: n_spec
  integer::niter

    call set_rates(t_rad_spec,h_rad_spec,aexp)

    ! Iteration to find mu
    err_mu=1.
    mu_left=0.5
    mu_right=1.3
    niter=0
    do while (err_mu > 1.d-4 .and. niter <= 50)
       mu_old=0.5*(mu_left+mu_right)
       !T = T2*mu_old
       call cmp_chem_eq(T,nH,t_rad_spec,n_spec,n_TOT,mu)
       err_mu = (mu-mu_old)/mu_old
       if(err_mu>0.)then
          mu_left =0.5*(mu_left+mu_right)
          mu_right=mu_right
       else
          mu_left =mu_left
          mu_right=0.5*(mu_left+mu_right)
       end if
       err_mu=ABS(err_mu)
       niter=niter+1
    end do
    if (niter > 50) then
       write(*,*) 'ERROR in calculation of mu : too many iterations.'
       STOP
    endif
end subroutine GetMuFromTemperature

