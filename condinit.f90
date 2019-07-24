!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use cooling_module, ONLY: twopi, kb, mh
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

  integer::i
  real(dp)::currad,xc,yc,zc,rho0g,rho0dm,mu,rho_cloud,c_s2,Phi0,PhiR,P_wind,P_cloud,nH
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_prs

  ! User-defined initial conditions

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
  scale_prs = scale_d * scale_v**2

  xc = x1_c*boxlen
  yc = x2_c*boxlen
  zc = x3_c*boxlen
  nH = n0g*scale_nH
  mu = 1.2
!  call GetMuFromTemperature(T_cloud,nH,mu)
  rho0g=n0g*mu
  rho0dm = gravity_params(1)
  Phi0 = -2.0*twopi*rho0dm*R_s**2
  c_s2 = kb*T_cloud/(mu*mh)/scale_v**2 !square of isothermal sound speed in cloud centre
  P_wind = ndens_wind*T_wind/scale_T2

  do i=1,nn
#if NDIM==3
    currad = dsqrt((x(i,1)-xc)**2 + (x(i,2)-yc)**2 + (x(i,3)-zc)**2)
#else
    currad = dsqrt((x(i,1)-xc)**2 + (x(i,2)-yc)**2)
#endif
    PhiR = Phi0*R_s*dlog(1+currad/R_s)/currad
    rho_cloud = rho0g*dexp(-(PhiR-Phi0)/c_s2)
    P_cloud = (rho_cloud*T_cloud/mu)/scale_T2
    if (P_cloud < P_wind) then
      q(i,1) = ndens_wind*mu
      q(i,ndim+2) = P_wind
    else
      q(i,1) = rho_cloud
      q(i,ndim+2) = P_cloud
    endif
    if (currad < Rad_cloud) then
      q(i,2) = 0.0      !x-velocity
      q(i,3) = 0.0      !y-velocity
#if NDIM==3
      q(i,4) = 0.0      !z-velocity
#endif
      q(i,imetal) = Z_cloud*0.02 !metallicity converted from relative to solar to absolute assuming Z_sol=0.02 as hardcoded in other parts of RAMSES
      q(i,imetal+1) = 1.0 !tracer
    else
      q(i,2) = 0.0      !x-velocity
      q(i,3) = vel_wind*1.e5/scale_v !y-velocity (given in km/s)
#if NDIM==3
      q(i,4) = 0.0      !z-velocity
#endif
      q(i,imetal) = Z_wind*0.02
      q(i,imetal+1) = 0.0	!tracer
    endif
!    write(*,*) "i ", i, " r ", currad, " rho: ", q(i,1), " P:",q(i,ndim+2), "mu:", mu,"Phi:",PhiR!" vy: ", q(i,3), " P: ", q(i,ndim+2), " T ", (q(i,ndim+2)*mu/q(i,1))*scale_t2, " x1: ", x1_c
  enddo
!write(*,*)'r_cut:',r_cut

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

