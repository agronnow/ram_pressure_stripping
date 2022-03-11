!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_commons, ONLY: t
  use amr_parameters
  use hydro_parameters
  use cooling_module, ONLY: kb
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
#else
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
#endif
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! If MHD, then:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E,
  ! U(i,6:8): Bleft, U(i,nvar+1:nvar+3): Bright
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar),save::q ! Primitive variables
  real(dp),save::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_prs
  real(dp)::Pwind,nH,vel,cosmo_time
  real(dp),save::vel0=-1d0
  real(dp),allocatable,save::tab_t(:),tab_vel(:)
  real(dp),save::dt
  integer::itab,ilun
  integer,save::ntab
  character(len=256)::fileloc
  real(dp),save::mu_wind
  logical::file_exists
  logical,save::firstcall=.true.
  logical,save::first_vel_max=.true.

!#ifdef SOLVERmhd
!  do ivar=1,nvar+3
!#else
!  do ivar=1,nvar
!#endif
!     do i=1,ncell
!        u(i,ivar)=boundary_var(ibound,ivar)
!     end do
!  end do

  if (firstcall) then
     call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
     nH = ndens_wind*scale_nH
     call GetMuFromTemperature(T_wind,nH,mu_wind)

     if (vel_wind > 0.0)then
        fileloc=trim(orbitfile)
        inquire(file=fileloc,exist=file_exists)
        ntab = 0
        if(file_exists) then
           open(newunit=ilun, file=fileloc)
           do
              read(ilun,*,end=10)
              ntab=ntab+1
           end do
10         rewind(ilun)
           allocate(tab_t(ntab))
           allocate(tab_vel(ntab))
           do i=1,ntab
              read(ilun,*)tab_t(i), tab_vel(i)
           end do
           close(ilun)
        endif
        dt = tab_t(2) - tab_t(1)
        cosmo_time = tinit_sim*3.154e16/scale_t-tbeg_wind
        itab = idint((cosmo_time-tab_t(1))/dt)+1 !Assume table is evenly spaced in t
        vel0 = (tab_vel(itab)*(tab_t(itab+1) - cosmo_time) + tab_vel(itab+1)*(cosmo_time - tab_t(itab)))/dt
     endif
     firstcall = .false.
  endif
  Pwind = ndens_wind*T_wind/scale_T2

  if (vel_wind > 0.0)then
     cosmo_time = t+tinit_sim*3.154e16/scale_t-tbeg_wind
     itab = idint((cosmo_time-tab_t(1))/dt)+1 !Assume table is evenly spaced in t
     vel = (tab_vel(itab)*(tab_t(itab+1) - cosmo_time) + tab_vel(itab+1)*(cosmo_time - tab_t(itab)))/dt
     !Boost injection velocity to get the correct velocity inside the volume at the front of the galaxy
     !velocity_multiplier must be calibrated for each potential
     vel = vel*(1d0+velocity_multiplier*((vel-vel0)/10.0))
     if(vel > vel_max)then
        if (first_vel_max)then
           write(*,*)'WARNING: velocity ',vel,' exceeds vel_max, set to ',vel_max
           first_vel_max=.false.
        endif
        vel=vel_max
     endif
!     write(*,*)"t, vel: ",cosmo_time,vel
  else !Static run
     vel = 0.0
  endif

  q(1:ncell,1) = ndens_wind*mu_wind	!density
  q(1:ncell,2) = 0.0	        !x-velocity
  q(1:ncell,3) = vel !vel_wind*1.e5/scale_v	!y-velocity (given in km/s)
#if NDIM == 3
  q(1:ncell,4) = 0.0    !z-velocity
#endif
  q(1:ncell,ndim+2) = Pwind	!pressure
  q(1:ncell,imetal) = Z_wind*0.02	!metallicity
  q(1:ncell,imetal+1) = 0.0	!tracer
  q(1:ncell,imetal+2) = 0.0	!sf tracer

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:ncell,1)=q(1:ncell,1)
  ! velocity -> momentum
  u(1:ncell,2)=q(1:ncell,1)*q(1:ncell,2)
#if NDIM>1
  u(1:ncell,3)=q(1:ncell,1)*q(1:ncell,3)
#endif
#if NDIM>2
  u(1:ncell,4)=q(1:ncell,1)*q(1:ncell,4)
#endif
  ! kinetic energy
  u(1:ncell,ndim+2)=0.0d0
  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+0.5*q(1:ncell,1)*q(1:ncell,2)**2
#if NDIM>1
  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+0.5*q(1:ncell,1)*q(1:ncell,3)**2
#endif
#if NDIM>2
  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+0.5*q(1:ncell,1)*q(1:ncell,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+q(1:ncell,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:ncell,ndim+2+ivar)=q(1:ncell,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+u(1:ncell,ndim+2+ivar)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     u(1:ncell,ivar)=q(1:ncell,1)*q(1:ncell,ivar)
  end do
#endif

end subroutine boundana
