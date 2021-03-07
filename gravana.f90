!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_commons, ONLY: t
  use amr_parameters
  use poisson_parameters
  use cooling_module, ONLY: twopi
  use incgamma
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
  real(dp)::r,rx,ry,rz,rho0,xmass,ymass,zmass,acc,r_max
#ifdef EINASTO
  logical, save::firstcall = .true.
  real(dp), save::gamma3n
#endif

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
     if (pot_grow_rate > 0.0)then
        r_max = min(max(r_cut, r_cut*(1d0+pot_grow_rate*(t-t_pot_grow_start))), r_tidal) ! Evolving r_cut
     else
        r_max = r_cut
     endif
     if (r < r_max) then
#ifdef EINASTO
       if (firstcall) then
         gamma3n = cmpgamma(3d0*ein_n)
         firstcall = .false.
       endif
       acc = -2d0*twopi*rho0*R_s**3*ein_n*(gamma3n - gammainc3n((r/R_s)**(1d0/(ein_n))))/r**2
#else
       ! NFW
       acc=-2d0*twopi*rho0*R_s**3*(dlog(1+r/R_s)-r/(r+R_s))/r**2
#endif
     else
       acc=0.0
     endif
     if (M_plummer > 0.0)acc = acc - M_plummer*r/(r_plummer**2 + r**2)**(3d0/2d0)
     f(i,1)=acc*rx/r
     f(i,2)=acc*ry/r
#if NDIM == 3
     f(i,3)=acc*rz/r
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
