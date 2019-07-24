!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters
  use cooling_module, ONLY: twopi
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::i
  real(dp)::r,rx,ry,rz,rho0,xmass,ymass,zmass,acc_nfw

  rho0 =gravity_params(1)
  xmass = x1_c*boxlen
  ymass = x2_c*boxlen
  zmass = x3_c*boxlen

  do i=1,ncell
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
#if NDIM == 3
     rz=x(i,3)-zmass
#else
     rz = 0.0
#endif
     r=sqrt(rx**2+ry**2+rz**2)
     if (r < r_cut) then
       acc_nfw=-2*twopi*rho0*R_s**3*(dlog(1+r/R_s)-r/(r+R_s))/r**2
     else
       acc_nfw=0.0
     endif
     f(i,1)=acc_nfw*rx/r
     f(i,2)=acc_nfw*ry/r
#if NDIM == 3
     f(i,3)=acc_nfw*rz/r
#endif
  end do

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
