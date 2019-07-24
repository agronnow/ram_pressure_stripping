!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,ndens_wind,T_wind,vel_wind,Z_wind
  use hydro_parameters, ONLY: nvar,nener,boundary_var,gamma,imetal
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
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_prs
  real(dp)::Pwind,mu

!#ifdef SOLVERmhd
!  do ivar=1,nvar+3
!#else
!  do ivar=1,nvar
!#endif
!     do i=1,ncell
!        u(i,ivar)=boundary_var(ibound,ivar)
!     end do
!  end do

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
  scale_prs = scale_d * scale_v**2

!  nWind = 1.d-4
!  Twind = 1.8d6
  mu = 1.34
  Pwind = ndens_wind*kb*T_wind/scale_prs

  q(1:,1) = ndens_wind*mu	!density
  q(1:ncell,2) = 0.0	        !x-velocity
  q(1:ncell,3) = vel_wind*1.e5/scale_v	!y-velocity (given in km/s)
#if NDIM == 3
  q(1:ncell,4) = 0.0    !z-velocity
#endif
  q(1:ncell,ndim+2) = Pwind	!pressure
  q(1:ncell,imetal) = Z_wind	!metallicity
  q(1:ncell,imetal+1) = 0.0	!tracer

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
