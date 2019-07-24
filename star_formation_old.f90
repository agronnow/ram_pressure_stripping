#if NDIM==3
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH , twopi
  use random
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
  integer,parameter::tag=1120
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled.
  ! It modifies hydrodynamic variables according to mass conservation
  ! and assumes an isothermal transformation...
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::t0,d0,d00,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star,ndebris_tot,ilun=10
  integer ::igrid,ix,iy,iz,ind,i,n,iskip,nx_loc,idim
  integer ::ntot,ntot_all,nstar_corrected,ncell
  logical ::ok_free
  real(dp)::d,x,y,z,u,v,w,e,prs,tg,zg
  real(dp)::mstar,dstar,tstar,nISM,nCOM,phi_t,phi_x,theta,sigs,scrit,b_turb,zeta
  real(dp)::T2,nH,T_poly,cs2,cs2_poly,trel,t_dyn,t_ff,tdec,uvar
  real(dp)::ul,ur,fl,fr,trgv,alpha0
  real(dp)::sigma2,sigma2_comp,sigma2_sole,lapld,flong,ftot,pcomp=0.3
  real(dp)::divv,divv2,curlv,curlva,curlvb,curlvc,curlv2
  real(dp)::birth_epoch,factG
  real(kind=8)::mlost_all,mtot_all
#ifndef WITHOUTMPI
  real(kind=8)::mlost,mtot
#endif
  real(kind=8)::PoissMean
  real(dp),parameter::pi=0.5*twopi
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min,d1,d2,d3,d4,d5,d6
  real(dp)::mdebris
  real(dp),dimension(1:nvector)::sfr_ff
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_cell2,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::ind_debris
  integer ,dimension(1:nvector,0:twondim)::ind_nbor
  logical ,dimension(1:nvector),save::ok,ok_new=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  character(LEN=256)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2,A,B,C,emag,beta,fbeta
#endif
#if NENER>0
  integer::irad
#endif
  real(dp)::dpass,Tpass,gpass,rho_sfr
  real(dp),save::tot_sf = 0.0
  real(dp)::ESN, rho_SN
  integer::nSN

  ! TODO: when f2008 is obligatory - remove this and replace erfc_pre_f08 below by
  ! the f2008 intrinsic erfc() function:
  real(dp) erfc_pre_f08

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return
  if(static)return

  if(verbose)write(*,*)' Entering SN feedback'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Star formation time scale from Gyr to code units
  ! SFR apply here for long lived stars only
  t0=t_star*(1d9*365.*24.*3600.)/scale_t
  trel=sf_trelax*1d6*(365.*24.*3600.)/scale_t

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
     nISM = MAX(nCOM,nISM)
  endif
  d0   = nISM/scale_nH
  d00  = n_star/scale_nH

  scale_m = scale_d*scale_l**ndim

  ! Initial star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star/mass_sph
  endif
  dstar=mstar/vol_loc

  factG = 1d0
  if(cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

  ! Birth epoch as proper time
  if(use_proper_time)then
     birth_epoch=texp
  else
     birth_epoch=t
  endif

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e-uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           prs = e/d
           !uold(ind_cell(i),5)=e/d
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do
     end do
  end do

! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  !------------------------------------------------
  ! Compute number of SNe in each cell
  !------------------------------------------------
  ntot=0
  ndebris_tot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
dpass = 0
Tpass = 0
gpass = 0
!tot_sf=0
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
           ! Temperature criterion
           do i=1,ngrid
              T2=prs*scale_T2*(gamma-1.0)
              nH=max(uold(ind_cell(i),1),smallr)*scale_nH
              T_poly=T2_star*(nH/nISM)**(g_star-1.0)
              T2=T2-T_poly
              if(T2>4e4) then
                ok(i)=.false.
              else
                Tpass=Tpass+1
              endif
           end do
        ! Geometrical criterion
        if(ivar_refine>0)then
           do i=1,ngrid
              d=uold(ind_cell(i),ivar_refine)
              if(d<=var_cut_refine) then
                ok(i)=.false.
              else
                gpass=gpass+1
              endif
           end do
        endif
        tot_sf = 0
        ! Calculate number of new stars in each cell using Poisson statistics
        !do i=1,ngrid
        !   nstar(i)=0
           ! Poisson mean
        !   rho_sfr = 10.0**(0.9+1.91*dlog10(uold(ind_cell(i),1)*0.02439303604))
        !   PoissMean=rho_sfr*scale_t/(3600*24*365.25)*vol_loc*dtnew(ilevel)*2d33/scale_m/mstar
        !   tot_sf = tot_sf + PoissMean
        !enddo

        do i=1,ngrid
           if(ok(i))then
              ! Compute mean number of events
              d=uold(ind_cell(i),1)
              mcell=d*vol_loc
              ! Free fall time of an homogeneous sphere
              tstar= .5427*sqrt(1.0/(factG*max(d,smallr)))
              if(.not.sf_virial) sfr_ff(i) = eps_star
              ! Gas mass to be converted into stars
              mgas=dtnew(ilevel)*(sfr_ff(i)/tstar)*mcell
              !d = d*0.02439303604 ! Convert cell density from amu cm^-3 to M_sol pc^-3
              ! Poisson mean
              rho_sfr = 10.0**(0.9+1.91*dlog10(d)) ! Volumetric SFR in M_sol yr^-1 kpc^-3 from Bacchini et al. 2019
              PoissMean=rho_SN*rho_sfr*scale_t/(3600*24*365.25)*vol_loc*dtnew(ilevel) ! Get expected number of SNe formed taking units into account
              !PoissMean=mgas/mstar
              !write(*,*) "PoissMean: ", PoissMean, " rho_sfr: ", rho_sfr, " mstar: ", mstar, " dt: ", dtnew(ilevel), " vol: ", vol_loc
              !write(*,*) "dt (yr): ", dtnew(ilevel)*scale_t/(3600*24*365.25), " vol(kpc^3): ", vol_loc*scale_l**ndim
              !if((trel>0.).and.(.not.cosmo)) PoissMean = PoissMean*min((t/trel),1.0)
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nSN)
              if (nSN > 0) then
                 do i=1,nSN
                    uold(ind_cell(i),5) = uold(ind_cell(i),5) + ESN/vol_loc
                    write(*,*) "SN explosion!"
                 end do
              endif
           endif
        enddo
     end do
  end do
!  write(*,*) "Tot sf: ", tot_sf
!write(*,*) "dens crit (d>", d0, "): ", dpass, " temp crit: ", Tpass, " refine crit: ", gpass


  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           !e=uold(ind_cell(i),5)*d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           !e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           !e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           !do irad=0,nener-1
           !   e=e+uold(ind_cell(i),inener+irad)
           !end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           !uold(ind_cell(i),5)=e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

end subroutine star_formation
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getnbor(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! the input cell.
  !-----------------------------------------------------------------
  integer::i,j,iok,ind
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok


  if(ilevel==1)then
     write(*,*) 'Warning: attempting to form stars on level 1 --> this is not allowed ...'
     return
  endif

  ! Get father cell
  do i=1,ncell
     ind_father(i,0)=ind_cell(i)
  end do

  ! Get father cell position in the grid
  do i=1,ncell
     pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
  end do

  ! Get father grid
  do i=1,ncell
     ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
  end do

  ! Get neighboring father grids
  call getnborgrids(ind_grid_father,igridn,ncell)

  ! Loop over position
  do ind=1,twotondim

     ! Select father cells that sit at position ind
     do j=0,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              igridn_ok(iok,j)=igridn(i,j)
           end if
        end do
     end do

     ! Get neighboring cells for selected cells
     if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

     ! Update neighboring father cells for selected cells
     do j=1,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              if(icelln_ok(iok,j)>0)then
                 ind_father(i,j)=icelln_ok(iok,j)
              else
                 ind_father(i,j)=ind_cell(i)
              end if
           end if
        end do
     end do

  end do


end subroutine getnbor
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function erfc_pre_f08(x)

! complementary error function
  use amr_commons, ONLY: dp
  implicit none
  real(dp) erfc_pre_f08
  real(dp) x, y
  real(kind=8) pv, ph
  real(kind=8) q0, q1, q2, q3, q4, q5, q6, q7
  real(kind=8) p0, p1, p2, p3, p4, p5, p6, p7
  parameter(pv= 1.26974899965115684d+01, ph= 6.10399733098688199d+00)
  parameter(p0= 2.96316885199227378d-01, p1= 1.81581125134637070d-01)
  parameter(p2= 6.81866451424939493d-02, p3= 1.56907543161966709d-02)
  parameter(p4= 2.21290116681517573d-03, p5= 1.91395813098742864d-04)
  parameter(p6= 9.71013284010551623d-06, p7= 1.66642447174307753d-07)
  parameter(q0= 6.12158644495538758d-02, q1= 5.50942780056002085d-01)
  parameter(q2= 1.53039662058770397d+00, q3= 2.99957952311300634d+00)
  parameter(q4= 4.95867777128246701d+00, q5= 7.41471251099335407d+00)
  parameter(q6= 1.04765104356545238d+01, q7= 1.48455557345597957d+01)

  y = x*x
  y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6) + p5/(y+q5)+p4/(y+q4)+p3/(y+q3) &
       &       + p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
  if (x < ph) y = y+2d0/(exp(pv*x)+1.0)
  erfc_pre_f08 = y

  return

end function erfc_pre_f08
