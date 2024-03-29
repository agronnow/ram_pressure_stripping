!------------------------------------------------------------------------------------------------------------
! Supernova feedback without star particles by Asger Gronnow
! SNII events are added randomly weighted by SFR in each cell calculated from the gas mass
! SNIa events are added randomly weighted by the SFH and stellar distribution of the galaxy
! SFH for cosmic times before the start of the simulation is read from a table,
! this is appended by the integrated SFR of the galaxy calculated during the simulation
! Radial stellar distribution is read from a table and is assumed to be unchanging.
! SN injection is done according to a scheme selected in the parameter file, choices are:
! - Direct thermal feedback scheme of Joung & Mac Low (2006) (requires high resolution)
!    Optionally delayed cooling can be applied when necessary
!    Injection region has a radius that encloses approximately a specific mass that
!    ensures high initial temperatures to avoid overcooling (but at least 3 cell radius).
!    Mass is evenly redistributed within the injection region
! - Kinetic feedback scheme largely similar to Hopkins et al. (2018) and Gentry, Madau & Krumholz (2020)
!   but without mass and metal injection. Injection region is the immediately adjacent cells (3x3x3).
!   Also works when adjacent cells are on different refinement levels
! - Kinetic feedback scheme of Simpson et al. (2015)
!   BEWARE: This scheme is largely untested and probably not currently correctly implemented!
!------------------------------------------------------------------------------------------------------------

!#if NDIM==3
subroutine subgrid_sn_feedback(ilevel, icount)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use sn_feedback_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH , twopi, kb
  use random
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
  integer,parameter::tag=1120
#endif
  integer::ilevel
  integer::icount
  ! local constants
  !real(dp)::t0,d0,d00,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,ilun=10
  integer ::igrid,ix,iy,iz,ind,i,n,iskip,nx_loc,idim
!  integer ::ntot,ntot_all,nstar_corrected,ncell
  logical ::ok_free
  real(dp)::d,x,y,z,u,v,w,e,prs,tg,zg
  !real(dp)::mstar,dstar,tstar,nISM,nCOM,phi_t,phi_x,theta,sigs,scrit,b_turb,zeta
  real(dp)::T2,nH!,T_poly,cs2,cs2_poly,trel,t_dyn,t_ff,tdec,uvar
!  real(dp)::ul,ur,fl,fr,trgv,alpha0
!  real(dp)::sigma2,sigma2_comp,sigma2_sole,lapld,flong,ftot,pcomp=0.3
!  real(dp)::divv,divv2,curlv,curlva,curlvb,curlvc,curlv2
!  real(dp)::birth_epoch,factG
!  real(kind=8)::mlost_all,mtot_all
!#ifndef WITHOUTMPI
!  real(kind=8)::mlost,mtot
!#endif
  real(kind=8)::PoissMean,maxPoissMean,maxPoissMean_all
  real(dp),parameter::pi=0.5*twopi
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min,d1,d2,d3,d4,d5,d6
!  real(dp)::mdebris
!  real(dp),dimension(1:nvector)::sfr_ff
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_cell2!,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
!  integer ,dimension(1:nvector),save::ind_debris
  integer ,dimension(1:nvector,0:twondim)::ind_nbor
  logical ,dimension(1:nvector),save::ok,ok_new=.true.
!  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  character(LEN=256)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2,A,B,C,emag,beta,fbeta
#endif
#if NENER>0
  integer::irad
#endif
  real(dp)::rho_sfr
  real(dp),save::tot_sf = 0.0
  real(dp)::ESN
  logical,save::nosn = .true.
  integer::Tpass,Tfail,Tpass_all,Tfail_all
  integer,parameter::RADCELL_MAX=1


#ifndef WITHOUTMPI
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all
  real(dp),dimension(:,:),allocatable::xSN_all
  integer,dimension(:),allocatable::levelSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles.
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip
  integer::nSN,nSN_loc,nSN_tot,iSN
  integer,dimension(1:ncpu)::nSN_icpu
!  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::current_time
!  real(dp)::dx_min,vol_min
!  integer::nx_loc
!  integer,dimension(:),allocatable::ind_grid
!  logical,dimension(:),allocatable::ok_free
  integer,dimension(:),allocatable::SNfinestlevel,levelSN,levelSN_loc,ncellsSN
  real(dp),dimension(:),allocatable::mSN,mSN_loc,rSN,volSN,wtot
  real(dp),dimension(:,:),allocatable::xSN,xSN_loc
  logical,dimension(:),allocatable::SNcooling

  logical,save::firstcall = .true.
  logical,save::calc_sfr = .false.
  real(dp),save::t_sfrlog = 0.0
  real(dp)::cursfr = 0.0
  real(dp)::weight = 0.0
#ifdef SNIA_FEEDBACK
#define NPDFBINS 1000
#define NSFHISTMAX 10000
  real(dp),dimension(1:MAXLEVEL),save::sfr_tot = 0.0
  real(dp),dimension(1:NPDFBINS),save::CDF_SNIa = 0.0
  real(dp),dimension(1:NPDFBINS)::PDF_SNIa = 0.0
  real(dp),dimension(1:NDIM),save::xCloud = 0.0
  real(dp),dimension(1:NSFHISTMAX),save::sfhist = 0.0
  real(dp),dimension(1:NSFHISTMAX),save::t_sfhist = 0.0
  real(dp),save::binwidth = 0.0
  character(len=255)::dummyline
  integer,save::nhist = 1
  logical,save::sfhist_update = .false.
  logical,save::pot_rad_max = .false.
  character(len=255)::sniadist_fname
  real(dp)::unif_rand,r2,rho_SNIa,area,DTD_A,DTD_s,currad,r,potrad
  real(dp)::PoissMeanIa,c_s2,nHc,ctime,csfh,diff,mindiff,dt,sfr,sfr_tot_level,rcosphi,sndist,theta
  integer::nSNIa,nt,imin,stat,clevel
  real(dp),dimension(:,:),allocatable::xpdf,xSNIa,min_r2,min_r2_all
  integer,dimension(:),allocatable::levelSNIa
  logical::doSNIa
#endif

#ifdef DELAYED_SN
  real(dp),dimension(NMAX_SN,3)::sn_coords_new
  real(dp),dimension(NMAX_SN)::rsn_sq_new
  integer,dimension(NMAX_SN)::sn_level_new
  integer,dimension(NMAX_SN)::sn_isrefined_all
  integer::nrem,inew
#endif
  integer,save::nSN_alltime = 0
  integer::pdevsn2

  ! TODO: when f2008 is obligatory - remove this and replace erfc_pre_f08 below by
  ! the f2008 intrinsic erfc() function:
  real(dp) erfc_pre_f08

!if (ilevel/=levelmin)return

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  !if(ndim.ne.3)return
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
  vol_loc=dx_loc**3
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**3

  ESN=1d51/(10.*2d33)/scale_v**2

  scale_m = scale_d*scale_l**3

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

!  if (myid==200)write(*,*)'seed fb:',localseedsn
  ! If necessary, initialize random number generator
  if(localseedsn(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseedsn=allseed(myid,1:IRandNumSize)
  end if

! get values of uold for density and velocities in virtual boundaries
!#ifndef WITHOUTMPI
!  do ivar=1,nvar
!     call make_virtual_fine_dp(uold(1,ivar),ilevel)
!  end do
!#endif

  allocate(xSN_loc(1:1000,1:ndim),mSN_loc(1:1000),levelSN_loc(1:1000))

#ifdef SNIA_FEEDBACK
#ifdef SNIA_HERINGER19
  DTD_A = 10.0**(-12.15) ! +0.10 dex -0.13 dex, [M_sun^-1 yr^-1], from Heringer et al. (2019)
  DTD_s = -1.34          ! +0.19 -0.17, from Heringer et al. (2019)
#else
  DTD_A = 2.11d-13       ! [M_sun^-1 yr^-1], from Maoz et al. (2012)
  DTD_s = -1.12          ! From Maoz et al. (2012)
#endif

  doSNIa = .false.
#if defined(SNIA_FEEDBACK) && !defined(SN_INJECT)
  if (ilevel == nlevelmax) then
     doSNIa = .true.
! Only generate SNIa at the last subcycle of the finest level
!     if (icount == 2) then
!        if(nsubcycle(ilevel)==2)doSNIa = .true.
!     else
!        doSNIa = .true.
!     endif
  endif
#endif

  nSNIa = 0
  if ((ilevel == levelmin) .and. (calc_sfr)) then
    calc_sfr = .false.
    if (sfhist_update) then
       nhist = nhist + 1
       t_sfhist(nhist) = (t-tbeg_wind)*scale_t/3.154e16 + tinit_sim ! Cosmic time in Gyr
       sfhist(nhist) = sum(sfr_tot)
#if NDIM==2
       sfhist(nhist) = sfhist(nhist)*4.0*Rad_cloud/3.0
#endif
       sfr_tot = 0.0
       sfhist_update = .false.
       if (myid==1) then
          fileloc=trim(output_dir)//'snIa_sfr.log'
          ilun=130
          inquire(file=fileloc,exist=file_exist)
          if(.not.file_exist) then
             open(ilun, file=fileloc, form='formatted')
             write(ilun,*)"Time                      sum(sfr_tot)"
          else
             open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
          endif
          write(ilun,'(3E26.16)') t, sfhist(nhist)
          close(ilun)

          fileloc=trim(output_dir)//'sfr.log'
          ilun=130
          inquire(file=fileloc,exist=file_exist)
          if(.not.file_exist) then
             open(ilun, file=fileloc, form='formatted')
             write(ilun,*)"Time                      sum(sfr_tot)"
          else
             open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
          endif
          write(ilun,'(3E26.16)') t, cursfr
          close(ilun)
       endif
    else
       t_sfrlog = t*scale_t/3.154e16 + tinit_sim
       cursfr = sum(sfr_tot)
#if NDIM==2
       cursfr = cursfr*4.0*Rad_cloud/3.0
#endif
       sfr_tot = 0.0
       if (myid==1) then
          fileloc=trim(output_dir)//'sfr.log'
          ilun=130
          inquire(file=fileloc,exist=file_exist)
          if(.not.file_exist) then
             open(ilun, file=fileloc, form='formatted')
             write(ilun,*)"Time                      sum(sfr_tot)"
          else
             open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
          endif
          write(ilun,'(3E26.16)') t, cursfr
          close(ilun)
       endif
    endif
  endif

  if (firstcall) then
     ! Construct initial Star Formation History
     ! SFH table is assumed to use cosmic age (t=0 is current time, t=13.88 Gyr is big bang) and is reordered to use cosmic time (t=0 is big bang, t=tinit_sim is t_sim=0)
     fileloc=trim(sfhistfile)
     ilun=150
     inquire(file=fileloc,exist=file_exist)
     nhist = 0
     if(file_exist) then
        if (vel_wind > 0.0)sfr_boost = 1.0
        open(ilun, file=fileloc)
        read(ilun,*)dummyline
        do
           read(ilun,*, iostat=stat)ctime,csfh
           if (stat /= 0)exit
           if (ctime < 13.88-tinit_sim)cycle
           t_sfhist(nhist+1) = 13.88 - ctime
           if (t_sfhist(nhist+1) < 0.1)t_sfhist(nhist+1) = 0.1
           sfhist(nhist+1) = csfh*sfr_boost
           nhist=nhist+1 
        end do
!        t_sfhfile = t_sfhist(nhist)
        close(ilun)
        !Reverse arrays to order according to cosmic time (t=0 is big bang) rather than age (t=0 is current time)
        t_sfhist = t_sfhist(nhist:1:-1)
        sfhist = sfhist(nhist:1:-1)
        if (myid==1)then
          do i=1,nhist
            write(*,*)"SFH: ",i,t_sfhist(i),sfhist(i)
          enddo
        endif
     endif

     ! Append the SFH computed during previous runs, relevant if restarting an earlier run
     fileloc=trim(output_dir)//'snIa_sfr.log'
     inquire(file=fileloc,exist=file_exist)
     if(file_exist) then
        open(ilun, file=fileloc)
        read(ilun,*)dummyline
        do
           read(ilun,*, iostat=stat)ctime,csfh
           ctime = ctime*scale_t/3.154e16 + tinit_sim
           if ((stat /= 0) .or. (ctime > tinit_sim + (t-tbeg_wind)*scale_t/3.154e16))exit
           t_sfhist(nhist+1) = ctime
           if (t_sfhist(nhist+1) < 0.1)t_sfhist(nhist+1) = 0.1
           sfhist(nhist+1) = csfh
           nhist=nhist+1
           if(myid==1)write(*,*)"restart SFH: ",nhist,t_sfhist(nhist),sfhist(nhist)
        end do
!        t_sfhfile = t_sfhist(nhist)
        close(ilun)
     endif
     
     xCloud(1) = x1_c
     xCloud(2) = x2_c
#if NDIM==3
     xCloud(3) = x3_c
#endif
  endif
  
  if ((myid==1) .and. firstcall)then !((pot_grow_rate > 0.0) .or. firstcall) .and. .not.(pot_rad_max)) then
    ! Calculate SNIa radius Comulative Distribution Function for the current potential truncation radius
    potrad = 2d0*r_cut !r_cut*(1d0+pot_grow_rate*(t-t_pot_grow_start))
    if (potrad > r_tidal)then
      potrad = r_tidal
      pot_rad_max = .true. ! Potential truncation radius will not grow further, stop recalculating SNIa CDF
      sniadist_fname = 'snIa_pdf_final.dat'
    else
      sniadist_fname = 'snIa_pdf_init.dat'
    endif
    area = 0.0
    binwidth = potrad/(NPDFBINS-1.0)
    
    do i=1,NPDFBINS
       currad = (i-1)*binwidth
       PDF_SNIa(i) = (3.0/(2.0*twopi*r_plummer**3))*(1.0+currad**2/r_plummer**2)**(-2.5)
       area = area + PDF_SNIa(i)*binwidth
    enddo
    do i=1,NPDFBINS
       ! Normalize PDF
       PDF_SNIa(i) = PDF_SNIa(i)/area
    enddo
    do i=1,NPDFBINS
       CDF_SNIa(i) = sum(PDF_SNIa(1:i))*binwidth
    enddo
!#ifdef DEBUG_SNIA
!    if ((t==0.0) .or. pot_rad_max)then
       fileloc=trim(output_dir)//sniadist_fname
       ilun=130
       open(ilun, file=fileloc, form='formatted')
       write(ilun,*)"rad                      PDF                   CDF"
       do i=1,NPDFBINS
          write(ilun,'(3E26.16)') (i-1)*binwidth, PDF_SNIa(i), CDF_SNIa(i)
       enddo
       close(ilun)
       write(*,*)"Wrote ",sniadist_fname
!    endif
!#endif
  endif

  if ((ilevel == levelmin) .and. .not.(firstcall)) then
    if ((t-tbeg_wind)*scale_t/3.154e16 + tinit_sim - t_sfhist(nhist) > dt_sfhist) then
       if (myid==1)write(*,*)'Update sfh',(t-tbeg_wind)*scale_t/3.154e16,tinit_sim,t_sfhist(nhist),dt_sfhist
       sfhist_update = .true.
       calc_sfr = .true.
    else if ((t-tbeg_wind)*scale_t/3.154e16 + tinit_sim - t_sfrlog > dt_sfrlog) then
       if(myid==1)write(*,*)'Calculate sfr'
       calc_sfr = .true.
    endif
  endif

  firstcall = .false.

  if (doSNIa) then
    ! Find number of SNIa and generate random positions according to PDF on root processor then broadcast results to all cpus
    if(myid==1)then
      ! Calculate global SNIa rate as SNR = int_0^t SFR(t-t')*DTD(t') dt' with DTD(t')=A*(t'/Gyr)^s (see e.g. Heringer et al. 2019)
      rho_SNIa = 0.0
      do nt=1,nhist
         if (nt < nhist) then
            dt = t_sfhist(nt+1) - t_sfhist(nt)
         else
            dt = (t-tbeg_wind)*scale_t/3.154e16 + tinit_sim - t_sfhist(nhist)
         endif
         rho_SNIa = rho_SNIa + sfhist(nhist-nt+1)*DTD_A*(t_sfhist(nt))**DTD_s*dt*1d9 ! SNIa per year
         !write(*,*)nt,nhist-nt+1,rho_SNIa,t_sfhist(nt),sfhist(nhist-nt+1),dt,DTD_A,DTD_s
      end do
      PoissMeanIa=rho_SNIa*scale_t/(3600*24*365.25)*dtnew(ilevel)/SN_batch_size ! Get expected number of SNIa formed taking units into account
#if NDIM==2
      PoissMeanIa = PoissMeanIa*4.0*Rad_cloud/3.0
#endif
#ifdef DEBUG_SNIA
      fileloc=trim(output_dir)//'snIa.log'
      ilun=130
      inquire(file=fileloc,exist=file_exist)
      if(.not.file_exist) then
         open(ilun, file=fileloc, form='formatted')
         write(ilun,*)"Time                      rho_SNIa                   PoissMeanIa"
      else
         open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
      endif
      write(ilun,'(3E26.16)') (t-tbeg_wind), rho_SNIa, PoissMeanIa
      close(ilun)
#endif
      ! Compute Poisson realisation
      call poissdev(localseedsn,PoissMeanIa,nSNIa)
      if (nSNIa > 0) then
         allocate(xpdf(1:nSNIa,1:ndim),xSNIa(1:nSNIa,1:ndim))
         ! Generate random positions for SNIa according to PDF
         do iSN=1,nSNIa
            ! Generate random real on [0,1] and convert to value drawn from SNIa distance PDF
            call ranf(localseedsn, unif_rand)
            r = unif_rand
            mindiff = 1.e22
            imin = 1
            do i=1,NPDFBINS
               diff = abs(CDF_SNIa(i) - r)
               if (diff < mindiff) then
                  mindiff = diff
                  imin = i
               endif
            enddo
            sndist = binwidth*(imin-1)
            ! Generate random coordinates at this distance by picking random x, y, z on the surface of a sphere
            ! This ONLY works in 3D
            call ranf(localseedsn, unif_rand)
            theta = twopi*unif_rand ! random theta in [0, 2pi]
            call ranf(localseedsn, unif_rand)
            rcosphi = sndist*(2d0*unif_rand - 1d0) ! random R*cos(phi) in [-R, R]
            xpdf(iSN,1) = sqrt(sndist**2 - rcosphi**2)*cos(theta) + xCloud(1)*boxlen
            xpdf(iSN,2) = sqrt(sndist**2 - rcosphi**2)*sin(theta) + xCloud(2)*boxlen
            xpdf(iSN,3) = rcosphi + xCloud(3)*boxlen
            write(*,*)'nSNIa',nSNIa,'SNIa',iSN,'rad',sndist,'theta',theta,'rcosphi',rcosphi,'r',r,'x,y,z=',xpdf(iSN,1),xpdf(iSN,2),xpdf(iSN,3)
         enddo
      endif
    endif

    call MPI_BCAST(nSNIa,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    if (nSNIa > 0) then
      if (myid/=1)allocate(xpdf(1:nSNIa,1:ndim),xSNIa(1:nSNIa,1:ndim))
      allocate(min_r2(1:2,1:nSNIa),levelSNIa(1:nSNIa))
      call MPI_BARRIER(MPI_COMM_WORLD,info)
      call MPI_BCAST(xpdf,nSNIa*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
      call MPI_BCAST(xSNIa,nSNIa*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
      call MPI_BCAST(PoissMeanIa,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)

      do iSN=1,nSNIa
        min_r2(1,iSN) = 1.e22 ! Huge value, must always be greater than min distance to SN
        min_r2(2,iSN) = myid  ! Tag with rank on each cpu for later use with MPI_MINLOC
      enddo
      if(myid==1)write(*,*) 'SNIa explosion'

     ! Find location of SNIa, i.e. the cell closest to the randomly generated position
     ! Loop over levels
     do clevel=levelmin,nlevelmax
        ! Cells center position relative to grid center position
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*0.5D0**clevel
           xc(ind,2)=(dble(iy)-0.5D0)*0.5D0**clevel
#if NDIM==3
           xc(ind,3)=(dble(iz)-0.5D0)*0.5D0**clevel
#endif
        end do

        ! Loop over grids
        ncache=active(clevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(clevel)%igrid(igrid+i-1)
           end do

           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax

              do i=1,ngrid
                 if(son(iskip+ind_grid(i))==0)then
                     do iSN=1,nSNIa
                         x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                         y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                         r2 = (x-xpdf(iSN,1))*(x-xpdf(iSN,1))+(y-xpdf(iSN,2))*(y-xpdf(iSN,2))
#if NDIM==3
                         z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                         r2 = r2 + (z-xpdf(iSN,3))*(z-xpdf(iSN,3))
#endif
                         if (r2 < min_r2(1,iSN)) then
                            min_r2(1,iSN)=r2
                            xSNIa(iSN,1) = x
                            xSNIa(iSN,2) = y
#if NDIM==3
                            xSNIa(iSN,3) = z
#endif
                            levelSNIa(iSN) = clevel
                         endif
                     enddo
                 endif
              enddo
           enddo
        enddo
      enddo
    endif
  endif
    
#endif

  !------------------------------------------------
  ! Compute number of SNe in each cell
  !------------------------------------------------

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  !dpass = 0
  Tpass = 0
  !gpass = 0
  Tfail=0
  maxPoissMean = 0.0
#ifdef SNIA_FEEDBACK
  sfr = 0.0
#endif

  pdevsn2=0
! get values of uold for density and velocities in virtual boundaries
!#ifndef WITHOUTMPI
!  do ivar=1,nvar
!     call make_virtual_fine_dp(uold(1,ivar),ilevel)
!  end do
!#endif

  nSN_loc = 0
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
           if (ok(i))then
              d=max(uold(ind_cell(i),1),smallr)
              u=uold(ind_cell(i),2)/d
              v=uold(ind_cell(i),3)/d
#if NDIM==3
              w=uold(ind_cell(i),4)/d
#else
              w=0.0
#endif
              e=uold(ind_cell(i),ndim+2)
              prs=e
#ifdef SOLVERmhd
              bx1=uold(ind_cell(i),6)
              by1=uold(ind_cell(i),7)
#if NDIM==3
              bz1=uold(ind_cell(i),8)
#endif
              bx2=uold(ind_cell(i),nvar+1)
              by2=uold(ind_cell(i),nvar+2)
#if NDIM==3
              bz2=uold(ind_cell(i),nvar+3)
              prs=prs-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#else
              prs=prs-0.125d0*((bx1+bx2)**2+(by1+by2)**2)
#endif
#endif
              prs=prs-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=0,nener-1
                 prs=prs-uold(ind_cell(i),inener+irad)
              end do
#endif
              prs = prs*(gamma-1.0)
              T2=prs*scale_T2*0.6/d
              !           nH=max(uold(ind_cell(i),1),smallr)*scale_nH
              !          T_poly=T2_star*(nH/nISM)**(g_star-1.0)
              !          T2=T2-T_poly
              if(T2>4e4) then
                 ok(i)=.false.
                 Tfail=Tfail+1
              else
                 Tpass=Tpass+1
              endif
           endif
        end do
!        ! Geometrical criterion
!        if(ivar_refine>0)then
!           do i=1,ngrid
!              d=uold(ind_cell(i),ivar_refine)
!              if(d<=var_cut_refine) then
!                ok(i)=.false.
!              else
!                gpass=gpass+1
!              endif
!           end do
!        endif
        tot_sf = 0

        !SNII
        do i=1,ngrid
           if(ok(i))then !Compute number of SNII if temperature is low enough for star formation
#ifndef SN_INJECT
              d=uold(ind_cell(i),1) 
              d = d*0.02439303604/1.36 ! Convert cell density from H+He density in amu cm^-3 to hydrogen density in M_sol pc^-3 assuming a helium fraction of 0.36 as in Gatto et al. 2013
              if (d > 0d0) then
                rho_sfr = vsfr_fac*d**vsfr_pow ! 10.0**(0.9+1.91*dlog10(d)) ! Volumetric SFR in M_sol yr^-1 kpc^-3 from Bacchini et al. 2019
#ifdef SNIA_FEEDBACK
                if (calc_sfr) then
                  sfr = sfr + rho_sfr*vol_loc
!                  write(*,*)'sfr',sfr,'rho_sfr',rho_sfr,'vol_loc',vol_loc,'d',d
                endif
#endif
                if(allow_coarse_sn .or. (ilevel==nlevelmax))then
                ! Poisson mean
                PoissMean=rho_SN*rho_sfr*scale_t/(3600*24*365.25)*vol_loc*dtnew(ilevel)/SN_batch_size ! Get expected number of SNe formed taking units into account
#if NDIM==2
                PoissMean = PoissMean*4.0*Rad_cloud/3.0
#endif
                ! Compute Poisson realisation
!                oldseed = localseedsn
                call poissdev(localseedsn,PoissMean,nSN)
                if (myid==200)pdevsn2=pdevsn2+1
                if (PoissMean > maxPoissMean) maxPoissMean = PoissMean
!               if (nosn) then
#else
!SN injection test sim
                if ((nSN_alltime==0).and.(t<1d-3))then
                   x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                   y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                   z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
!                   write(*,*)abs(x - 0.51*boxlen), abs(y - 0.5*boxlen), abs(z - 0.5*boxlen), dx_loc
                   if ((abs(x - SN_inject_x) < 0.5*dx_min) .and. (abs(y - SN_inject_y) < 0.5*dx_min) .and. (abs(z - SN_inject_z) < 0.5*dx_min))then
!if ((abs(x - 0.5*boxlen < 0.51*dx_min) .and. (x > 0.5*boxlen) .and. (y - 0.5*boxlen < 0.51*dx_min) .and. (z - 0.5*boxlen < 0.51*dx_min) .and. (y > 0.5*boxlen) .and. (z > 0.5*boxlen))then
                      nSN = 1
                      nSN_alltime = 1
                   else
                      nSN = 0
                   endif
                else
                   nSN = 0
                endif
                PoissMean = 0.0
                rho_sfr = 0.0
#endif
                do iSN=1,nSN
                   ! Get gas cell position
                   x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                   y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
#if NDIM==3
                   z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
#endif
                   if (y > x2_c*boxlen + 3.0)then !Skip SNe in stripped gas >3 kpc behind galaxy
                      write(*,*)"Skipping SN far upstream at ",x,y,z
                   else
                      nSN_loc=nSN_loc+1
                      xSN_loc(nSN_loc,1)=x
                      xSN_loc(nSN_loc,2)=y
#if NDIM==3
                      xSN_loc(nSN_loc,3)=z
#endif
                      levelSN_loc(nSN_loc) = ilevel
                      mSN_loc(nSN_loc)=SN_batch_size*10.0*2d33/(scale_d*scale_l**3) !Always assume 10 solar mass ejection

                      fileloc=trim(output_dir)//'sn.dat'
                      ilun=140
                      inquire(file=fileloc,exist=file_exist)
                      if(.not.file_exist) then
                         open(ilun, file=fileloc, form='formatted')
                         write(ilun,*)"Time                      x                         y                         z                          Level ProbSN                    rho_sfgas                 ProcID"!  Seeds"
                      else
                         open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
                      endif
#if NDIM==2
                      z = 0.0
#endif
!#ifdef SN_INJECT
                      write(ilun,'(4E26.16,I5,E26.16,E26.16,I5)') t, x, y, z, ilevel, PoissMean, uold(ind_cell(i),1), myid!, oldseed(1), oldseed(2), oldseed(3), oldseed(4) !No passive tracer in SN injection test sim
!#else
!                      write(ilun,'(4E26.16,I5,E26.16,E26.16,I5)') t, x, y, z, ilevel, PoissMean, uold(ind_cell(i),1), myid!, oldseed(1), oldseed(2), oldseed(3), oldseed(4)
!#endif
                      close(ilun)
                   endif
                !   do iSN=1,nSN
                !      uold(ind_cell(i),ndim+2) = uold(ind_cell(i),ndim+2) + 10.0*2d33/(scale_d*scale_l**3)*ESN/vol_loc
                !      write(*,*) "SN explosion!"
                !   end do
                end do
                endif
#ifndef SN_INJECT
              else
                rho_sfr = 0.0
                nSN = 0
              endif
#endif
           endif
        enddo
     end do
  end do

#ifdef SNIA_FEEDBACK
  if (nSNIa > 0) then
     !write(*,*)'min_r2',min_r2(1,1),xSNIa(1,1),xSNIa(1,2),xSNIa(1,3),myid
     allocate(min_r2_all(1:2,1:nSNIa))
     call MPI_ALLREDUCE(min_r2,min_r2_all,nSNIa,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,info)

     do iSN=1,nSNIa
        if (int(min_r2_all(2,iSN)) == myid) then ! The cell containing the SN center is on this processor
           if (levelSNIa(iSN) == nlevelmax)then
              nSN_loc=nSN_loc+1
              xSN_loc(nSN_loc,1)=xSNIa(iSN,1)
              xSN_loc(nSN_loc,2)=xSNIa(iSN,2)
#if NDIM==3
              xSN_loc(nSN_loc,3)=xSNIa(iSN,3)
#endif
              mSN_loc(nSN_loc)=SN_batch_size*10.0*2d33/(scale_d*scale_l**3) !Always assume 10 solar mass ejection
              levelSN_loc(nSN_loc) = levelSNIa(iSN)
              !if (outputSN) then
              fileloc=trim(output_dir)//'snIa.dat'
              ilun=140
              inquire(file=fileloc,exist=file_exist)
              if(.not.file_exist) then
                 open(ilun, file=fileloc, form='formatted')
                 write(ilun,*)"Time                      x                         y                         z ProbSN                    Level         ProcID"!  Seeds"
              else
                 open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
              endif
#if NDIM==3
              z = xSNIa(iSN,3)
#else
              z = 0.0
#endif
              write(ilun,'(5E26.16,E26.16,2I5)') t, xSNIa(iSN,1), xSNIa(iSN,2), z, sqrt((xSNIa(iSN,1)-x1_c*boxlen)**2+(xSNIa(iSN,2)-x2_c*boxlen)**2+(z-x3_c*boxlen)**2), PoissMeanIa, levelSNIa(iSN),myid!, oldseed(1), oldseed(2), oldseed(3), oldseed(4)
              close(ilun)
           else
              write(*,*)"WARNING: Skipping SNIa on coarse level: Level ",levelSNIa(iSN)," x,y,z: ",xSNIa(iSN,1), xSNIa(iSN,2), xSNIa(iSN,3)
           endif
        endif
     end do
  endif
#endif

  nSN_icpu=0
  nSN_icpu(myid)=nSN_loc
  Tpass_all = 0
  Tfail_all = 0
  maxPoissMean_all = 0.0

  sfr_tot_level = 0.0
#ifndef WITHOUTMPI
  if (calc_sfr) then
    call MPI_ALLREDUCE(sfr,sfr_tot_level,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,info)
    if(myid==1)write(*,*)'reduce sfr: sfr on cpu 1=',sfr,' sfr_tot=',sfr_tot_level,' level=',ilevel
  endif

  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all

  if (prob_debug) then
    call MPI_REDUCE(Tpass,Tpass_all,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,info)
    call MPI_REDUCE(Tfail,Tfail_all,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,info)
    call MPI_REDUCE(maxPoissMean,maxPoissMean_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,info)
    if (myid==1)write(*,*)'Tpass:',Tpass_all,' Tfail:', Tfail_all,' maxPoissMean:',maxPoissMean_all
  endif
#else
  sfr_tot_level = sfr
#endif

  if (calc_sfr) then
     weight = nsubcycle(levelmin)/(1d0*product(nsubcycle(levelmin:ilevel))) !weight by number of subcycles on level such that the SFR on finer levels are correctly averaged when finding the total SFR during next coarse time step
     sfr_tot(ilevel-levelmin+1) = sfr_tot(ilevel-levelmin+1) + weight*sfr_tot_level
  endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  !if(myid==200)write(*,*)'calls to poissdev: ',pdevsn2

#ifdef DELAYED_SN
  if ((nSN_tot .eq. 0).and.(nSN_prev .eq. 0)) return

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN explosions (possibly including ones that will be delayed)=',nSN_tot
     write(*,*)'Number of delayed SN explosions=',nSN_prev
     write(*,*)'-----------------------------------------------'
  endif
#else
  if (nSN_tot .eq. 0) return

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN explosions=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif
#endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN(1:nSN_tot,1:ndim))
  allocate(mSN(1:nSN_tot),levelSN(1:nSN_tot))
  xSN=0.;mSN=0.;levelSN=0
  ! Allocate arrays for particles index and parent grid
!  if(nSN_loc>0)then
!     allocate(ind_part(1:nSN_loc),ind_grid(1:nSN_loc),ok_free(1:nSN_loc))
!  endif

  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif

  ip=0
  do icpu=1,ncpu
    do i=1,nSN_loc
      xSN(i+iSN,1) = xSN_loc(i,1)
      xSN(i+iSN,2) = xSN_loc(i,2)
#if NDIM==3
      xSN(i+iSN,3) = xSN_loc(i,3)
#endif
     mSN(i+iSN) = mSN_loc(i)
     levelSN(i+iSN) = levelSN_loc(i)
    end do
  end do
  deallocate(xSN_loc,mSN_loc,levelSN_loc)

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:ndim),mSN_all(1:nSN_tot),levelSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN,xSN_all,nSN_tot*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(levelSN,levelSN_all,nSN_tot  ,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  xSN=xSN_all
  mSN=mSN_all
  levelSN=levelSN_all
  deallocate(xSN_all,mSN_all,levelSN_all)
#endif

  nSN=nSN_tot
  nSN_alltime = nSN_alltime + nSN
  
#ifdef DELAYED_SN
  if(nSN_prev > 0)then
     call MPI_ALLREDUCE(sn_isrefined,sn_isrefined_all,nSN_prev  ,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     sn_isrefined = sn_isrefined_all
     allocate(rSN(1:nSN_prev),SNfinestlevel(1:nSN_prev),volSN(1:nSN_prev),wtot(1:nSN_prev),ncellsSN(1:nSN_prev),SNcooling(1:nSN_prev))
     ! Add SN from previous time step
     if(myid==1)write(*,*)nSN_prev,'delayed SNe at ',sn_coords(1,1),' ',sn_coords(1,2),' ',sn_coords(1,3)
     ! Compute blast radius
     call subgrid_average_SN(sn_coords(1:nSN_prev,1:ndim),mSN,rSN,volSN,levelSN,wtot,ncellsSN,nSN_prev,SNfinestlevel,SNcooling,.true.)

     ! Modify hydro quantities to account for a Sedov blast wave
     call subgrid_Sedov_blast(sn_coords(1:nSN_prev,1:ndim),mSN,rSN,volSN,levelSN,wtot,ncellsSN,nSN_prev,SNfinestlevel,SNcooling,.true.)
     deallocate(rSN,SNfinestlevel,volSN,levelSN,wtot,ncellsSN,SNcooling)
  endif
#endif

  allocate(rSN(1:nSN),SNfinestlevel(1:nSN),volSN(1:nSN),wtot(1:nSN),ncellsSN(1:nSN),SNcooling(1:nSN))

  ! Compute blast radius
  call subgrid_average_SN(xSN,mSN,rSN,volSN,levelSN,wtot,ncellsSN,nSN,SNfinestlevel,SNcooling,.false.)

  ! Modify hydro quantities to account for a Sedov blast wave
  call subgrid_Sedov_blast(xSN,mSN,rSN,volSN,levelSN,wtot,ncellsSN,nSN,SNfinestlevel,SNcooling,.false.)

#ifdef DELAYED_SN
  inew=0
  nrem = 0
  if (nSN_prev > 0) then
     if (nSN_prev > NMAX_SN)write(*,*)"WARNING: Too many SN for DELAYED_SN module!!!"
     do iSN=1,nSN_prev
        if (sn_isrefined(iSN)==0)then
           inew = inew+1
           sn_level_new(inew) = sn_level(iSN)
           sn_coords_new(inew,1:3) = sn_coords(iSN,1:3)
           rsn_sq_new(inew) = rsn_sq(iSN)
        else
           nrem = nrem+1
        endif
     enddo
     nSN_prev = inew
     if (inew > 0)then
        sn_level(1:inew) = sn_level_new(1:inew)
        sn_coords(1:inew,1:3) = sn_coords_new(1:inew,1:3)
        rsn_sq(1:inew) = rsn_sq_new(1:inew)
        sn_isrefined(1:inew) = 0
        if (myid==1)then
           do iSN=1,inew
              write(*,*)"delayed SN: ",iSN," ",sn_coords(iSN,1)," ",sn_coords(iSN,2)," ",sn_coords(iSN,3)," ",sn_level(iSN)," ",rsn_sq(iSN)
           enddo
        endif
     endif
  endif
  if(myid==1)write(*,*)"Removed ",nrem," delayed SNe"
  do iSN=1,nSN
      if (SNfinestlevel(iSN) > 0) then
        nSN_prev = nSN_prev + 1
        sn_level(nSN_prev) = SNfinestlevel(iSN)
        sn_coords(nSN_prev,1:3) = xSN(iSN,1:3)
        rsn_sq(nSN_prev) = (2d0*rSN(iSN))**2
        sn_isrefined(nSN_prev) = 0
        if(myid==1)write(*,*)"nSNdelay",nSN_prev," SNcoords:",sn_coords(nSN_prev,1),",",sn_coords(nSN_prev,2),",",sn_coords(nSN_prev,3)
     endif
  enddo
#endif

  deallocate(xSN,mSN,rSN,volSN,wtot,levelSN,SNcooling)

#ifdef SNIA_FEEDBACK
  if (nSNIa > 0)deallocate(xpdf,xSNIa,min_r2,min_r2_all,levelSNIa)
#endif


  ! Update hydro quantities for split cells
!  do ilevel=nlevelmax,levelmin,-1
!     call upload_fine(ilevel)
!     do ivar=1,nvar
!        call make_virtual_fine_dp(uold(1,ivar),ilevel)
!     enddo
!  enddo

  if(verbose)write(*,*)'Exiting subgrid_sn_feedback'

end subroutine subgrid_sn_feedback
!#endif

!################################################################
!################################################################
!################################################################
!################################################################
subroutine subgrid_average_SN(xSN,mSN,rSN,SNvol,level_SN,wtot,ncellsSN,nSN,SNfinestlevel,SNcooling,delayed)
  use pm_commons
  use amr_commons
  use hydro_commons
#if NDIM==2
  use cooling_module, only: twopi
#endif
  use mpi_mod
  use sn_feedback_commons
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer,parameter::RADCELL_MAX=1
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,iSN,ind,ix,iy,iz,ngrid,iskip,radcells
  integer::i,nx_loc,igrid,ivar
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,u,v,w,u2,v2,w2,dr_cell,massdiff,mindiff,dprev,momprev,momnew,fZ,mcenter
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax,dx_SN,adjacency,cellweight
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc
  integer ,dimension(1:nSN)::SNfinestlevel,flagrefine,flagrefine_all,level_SN,ncellsSN
#ifndef DELAYED_SN
  integer ,dimension(1:nSN)::SNmaxrad,SNmaxrad_all
#endif
  real(dp),dimension(1:nSN)::ekBlast,rSN,vol_center,vol_center_all,wtot,wtot_all,mSN
  logical,dimension(1:nSN)::SNcooling
  real(dp),dimension(1:nSN,1:RADCELL_MAX)::vol_gas,vol_gas_all,mtot,mtot_all
  integer,dimension(1:nSN,1:RADCELL_MAX)::snmaxlevel,snmaxlevel_all,snncells,snncells_all,SNcoarsestlevel,SNcoarsestlevel_all
  real(dp),dimension(1:nSN,1:ndim)::xSN
  logical::file_exist,update_boundary
  integer::ilun
  character(LEN=256)::fileloc
  real(dp),dimension(1:nSN)::SNmenc,SNvol
  logical ,dimension(1:nvector),save::ok

  logical::delayed ! Don't flag SN for refinement and delayed blast multiple times
  logical::skip,kinetic_inj
  character(len=5)::delayedstr1,delayedstr,coolstr

  if(nSN==0)return
  if(verbose)write(*,*)'Entering average_SN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Initialize the averaged variables
  vol_gas=0.0;rSN=0.0;vol_gas_all=0.0;mtot=0.0;mtot_all=0.0;SNmenc=0.0;SNvol=0.0;ncellsSN=0
  snmaxlevel=0;snmaxlevel_all=0;flagrefine=0;flagrefine_all=0;snncells=0;snncells_all=0;wtot=0.0;wtot_all=0.0;SNcoarsestlevel=nlevelmax;SNcoarsestlevel_all=0
#ifndef DELAYED_SN
  SNmaxrad=RADCELL_MAX
#endif

  update_boundary = .false.

  do iSN=1,nSN
#ifdef DELAYED_SN
     if (delayed)then
        if(sn_isrefined(iSN)==0)cycle
     endif
#endif
        dx_SN = 0.5D0**level_SN(iSN)*scale
        ! Loop over levels
        do ilevel=levelmin,nlevelmax
           ! Computing local volume (important for averaging hydro quantities)
           dx=0.5D0**ilevel
           dx_loc=dx*scale
           vol_loc=dx_loc**3
#if NDIM==2
           vol_loc=2.0*twopi*vol_loc/3.0
#endif
           ! Cells center position relative to grid center position
           do ind=1,twotondim
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              xc(ind,1)=(dble(ix)-0.5D0)*dx
              xc(ind,2)=(dble(iy)-0.5D0)*dx
#if NDIM==3
              xc(ind,3)=(dble(iz)-0.5D0)*dx
#endif
           end do

           ! Loop over grids
           ncache=active(ilevel)%ngrid
           do igrid=1,ncache,nvector
              ngrid=MIN(nvector,ncache-igrid+1)
              do i=1,ngrid
                 ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
              end do

              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 do i=1,ngrid
                    ind_cell(i)=iskip+ind_grid(i)
                 end do

                 ! Flag leaf cells
                 do i=1,ngrid
                    ok(i)=son(ind_cell(i))==0
                 end do

                 do i=1,ngrid
                    if(ok(i))then
                       ! Get gas cell position
                       x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                       y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
#if NDIM==3
                       z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
#endif
                       ! Check if the cell lies within the SN radius
                       dxx=x-xSN(iSN,1)
                       dyy=y-xSN(iSN,2)
#if NDIM==3
                       dzz=z-xSN(iSN,3)
                       dr_SN=dxx**2+dyy**2+dzz**2
                       dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
#else
                       dr_SN=dxx**2+dyy**2
                       dr_cell=MAX(ABS(dxx),ABS(dyy))
#endif
                       do radcells=1,RADCELL_MAX
                          if(((dr_SN .lt. (dx_SN*(radcells+0.5))**2) .and. (radcells>1)) .or. ((dr_cell < 1d-9 + dx_SN/2d0 + dx_loc/2d0) .and. (min(abs(dxx),abs(dyy),abs(dzz)) < 1d-9))) then
!                             if ((ilevel ~= plevel) .and. (plevel >= 0) .and. (radcells <= SNmaxrad(iSN))) then
!                                 SNmaxrad(iSN) = radcells-1 ! SN radius must be smaller than this to avoid overlapping coarse cells
!                             endif
                             if (.not. delayed) then
                                if (ilevel > snmaxlevel(iSN,radcells))then
                                   if (snmaxlevel(iSN,radcells) > 0)flagrefine(iSN) = 1
                                   snmaxlevel(iSN,radcells) = ilevel
                                endif
                                if (ilevel < sncoarsestlevel(iSN,radcells))sncoarsestlevel(iSN,radcells) = ilevel
                             endif
                             mtot(iSN,radcells) = mtot(iSN,radcells) + max(uold(ind_cell(i),1),smallr)*vol_loc
                             vol_gas(iSN,radcells) = vol_gas(iSN,radcells) + vol_loc
                             snncells(iSN,radcells) = snncells(iSN,radcells) + 1
                             if(dr_SN < 1d-10)mcenter=max(uold(ind_cell(i),1),smallr)*vol_loc
!write(*,*)'radcells',radcells,' mtot',mtot(iSN,radcells),' vol',vol_gas(iSN,radcells)
                          endif
                       enddo
                    endif
                 end do
              end do
           end do     ! End loop over grids
        end do    ! End loop over levels
!        write(*,*)"mtot",mtot," rad",sqrt(rSN(iSN))
!        if (mtot < SN_blast_mass)rSN(iSN) = rSN(iSN) + dx_loc
  enddo

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas,vol_gas_all,nSN*RADCELL_MAX  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,nSN*RADCELL_MAX  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(snncells,snncells_all,nSN*RADCELL_MAX  ,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sncoarsestlevel,sncoarsestlevel_all,nSN*RADCELL_MAX  ,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  if (.not. delayed)then
     call MPI_ALLREDUCE(snmaxlevel,snmaxlevel_all,nSN*RADCELL_MAX  ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(flagrefine,flagrefine_all,nSN  ,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  endif
#endif

#ifndef DELAYED_SN
  do iSN=1,nSN
     do radcells=1,RADCELL_MAX
        if ((sncoarsestlevel_all(iSN,radcells) < level_SN(iSN)) .and. (SNmaxrad(iSN) > radcells))SNmaxrad(iSN) = radcells-1
     enddo
  enddo
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(SNmaxrad,SNmaxrad_all,nSN  ,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#endif

  do iSN=1,nSN
     skip=.false.
#ifdef DELAYED_SN
     if (delayed)then
        if(sn_isrefined(iSN)==0)cycle
     endif
#endif
     dx_SN = 0.5D0**level_SN(iSN)*scale
     mindiff = 1d10
     do radcells=1,RADCELL_MAX
        massdiff = abs(SN_blast_mass - mtot_all(iSN,radcells)*scale_d*scale_l**3/2d33)
!write(*,*)'radcells',radcells,' massdiff',massdiff,' mtot',mtot_all(iSN,radcells)*scale_d*scale_l**3/2d33,' vol_gas',vol_gas_all(iSN,radcells)
        if (massdiff < mindiff) then
           mindiff = massdiff
           rSN(iSN) = radcells*dx_SN
           if (.not. delayed) then
              SNfinestlevel(iSN) = maxval(snmaxlevel_all(iSN,1:radcells))
              if (flagrefine_all(iSN) == 0)SNfinestlevel(iSN) = 0 ! Disable SN based refinement if entire SN blast is on the same level
           endif
           SNmenc(iSN) = mtot_all(iSN,radcells)
           SNvol(iSN) = vol_gas_all(iSN,radcells)
           ncellsSN(iSN) = snncells_all(iSN,radcells)
        endif
        !write(*,*)"Ncells: ",radcells, " ", snncells_all(iSN,radcells)
     enddo
     if (SNmenc(iSN)*scale_d*scale_l**3/2d33 < 1d0)then
        if(myid==1)write(*,*)"SN enclosed mass less than 1 M_sun, skipping this SN!!!"
        rSN(iSN)=-1d0
        cycle
      endif
!     if (rSN(iSN) > 1.1*RADCELL_MAX*dx_SN)then
!        if (myid==1)write(*,*)"WARNING: Skipping SN with radius greater than the maximum number of cells allowed:", RADCELL_MAX
!        skip=.true.
!     endif
     if (delayed)SNfinestlevel(iSN) = 0
     if ((momentum_fb).and.(allow_coarse_SN))SNmaxrad_all(iSN)=RADCELL_MAX
#ifndef DELAYED_SN
     if (SNmaxrad_all(iSN)*dx_SN < rSN(iSN)) then
        if (myid==1)write(*,*)"WARNING: SN should extend ",rSN(iSN), " kpc but can only extend ",SNmaxrad_all(iSN)*dx_SN," kpc without overlapping coarser cells!!!"
        rSN(iSN) = SNmaxrad_all(iSN)*dx_SN
        SNmenc(iSN) = mtot_all(iSN,SNmaxrad_all(iSN))
        SNvol(iSN) = vol_gas_all(iSN,SNmaxrad_all(iSN))
        ncellsSN(iSN) = snncells_all(iSN,SNmaxrad_all(iSN))
     endif
     if(SNmaxrad_all(iSN) == 0)then
        rSN(iSN) = 0
        SNmenc(iSN) = mcenter
        SNvol(iSN) = dx_SN**3
        ncellsSN(iSN) = 1
     endif
#endif
    if ((rSN(iSN) == 0) .or. (SNmaxrad_all(iSN)==0))then
       if (myid==1)write(*,*)"WARNING: SN has radius 0!"
       skip=.true.
    endif

    if (((momentum_fb.and.(rSN(iSN) < 3d0*dx_SN)) .or. ((rSN(iSN) < 3d0*dx_SN).and.delayed_cooling)) .and. (.not. skip))then
       if (delayed_cooling)then
         SNcooling(iSN) = .false.
         kinetic_inj = .false.
         radcells = min(SNmaxrad_all(iSN), 3)
       else
         SNcooling(iSN) = .true.
         kinetic_inj = .true.
         if((myid==1).and.(SNmaxrad_all(iSN)<mominj_rad))write(*,*)"WARNING: SN should extend ",mominj_rad, " cells but has been shrunk to only extend ",SNmaxrad_all(iSN)," cells to avoid overlapping coarser cells!"
         radcells = min(SNmaxrad_all(iSN), mominj_rad)
       endif
       rSN(iSN) = radcells*dx_SN
       SNmenc(iSN) = mtot_all(iSN,radcells)
       SNvol(iSN) = dx_SN**3*19!vol_gas_all(iSN,radcells)
       ncellsSN(iSN) = snncells_all(iSN,radcells)
    else
       SNcooling(iSN) = .true.
       kinetic_inj = .false.
    endif

     if(myid==1) then
       fileloc=trim(output_dir)//'snblast.dat'
       ilun=140
       inquire(file=fileloc,exist=file_exist)
       if(.not.file_exist) then
          open(ilun, file=fileloc, form='formatted')
#ifdef DELAYED_SN
          write(ilun,*)"Time                      x                         y                         z                         Ejecta mass (Msun)        Blast radius (kpc)        Blast radius (cells)        Blast cells              Cooling length (kpc)       Blast mass (Msun)          Blast temperature (K)   Delayed  To be delayed"
#else
          write(ilun,*)"Time                      x                         y                         z                         Ejecta mass (Msun)        Blast radius (kpc)        Blast radius (cells)        Blast cells              Cooling length (kpc)       Blast mass (Msun)          Blast temperature (K)   Max rad (cells)  Cooling"
#endif
       else
          open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
       endif
#if NDIM==3
       z = xSN(iSN,3)
#else
       z = 0.0
#endif
       if (SNcooling(iSN))then
          coolstr = 'Y'
       else
          coolstr = 'N'
       endif
       fZ = 2d0
       if (Z_cloud < 0.01)fZ=Z_cloud**(-0.14)
#ifdef DELAYED_SN
       if (delayed)then
          delayedstr1 = 'Y'
       else
          delayedstr1 = 'N'
       endif
       if (SNfinestlevel(iSN)==0)then
          delayedstr = 'N'
       else
          delayedstr = 'Y'
       endif
       if (skip)delayedstr = 'SKIP'
       write(ilun,'(6E26.16,2I5,3E26.16,3A7)') t, xSN(iSN,1), xSN(iSN,2), z, mSN(iSN)/(2d33/(scale_d*scale_l**3)), rSN(iSN), int(rSN(iSN)/dx_SN), ncellsSN(iSN), 0.0284*(mSN(iSN)/(10d0*2d33/(scale_d*scale_l**3)))**(2d0/7d0)*(SNmenc(iSN)/SNvol(iSN))**(-3d0/7d0)*fZ, SNmenc(iSN)*scale_d*scale_l**3/2d33, (SN_batch_size*1d51*(gamma-1d0)/(SNmenc(iSN)*scale_d*scale_l**3))*(0.6*1.66e-24/1.3806e-16), coolstr, delayedstr1,delayedstr
#else
       if (skip)coolstr = 'SKIP'
       write(ilun,'(6E26.16,2I5,3E26.16,I5,A7)') t, xSN(iSN,1), xSN(iSN,2), z, mSN(iSN)/(2d33/(scale_d*scale_l**3)), rSN(iSN), int(rSN(iSN)/dx_SN), ncellsSN(iSN), 0.0284*(mSN(iSN)/(10d0*2d33/(scale_d*scale_l**3)))**(2d0/7d0)*(SNmenc(iSN)/SNvol(iSN))**(-3d0/7d0)*fZ, SNmenc(iSN)*scale_d*scale_l**3/2d33, (SN_batch_size*1d51*(gamma-1d0)/(SNmenc(iSN)*scale_d*scale_l**3))*(0.6*1.66e-24/1.3806e-16),SNmaxrad_all(iSN), coolstr
#endif
       close(ilun)
     endif
     if((skip).or.(kinetic_inj))cycle

#ifdef DELAYED_SN
     if (SNfinestlevel(iSN) > 0)cycle
#endif

     ! Evenly redistribute mass within SN injection region for fully thermal injection
     ! Loop over levels
     do ilevel=levelmin,nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx=0.5D0**ilevel
        dx_loc=dx*scale
        vol_loc=dx_loc**3
#if NDIM==2
        vol_loc=2.0*twopi*vol_loc/3.0
#endif
        ! Cells center position relative to grid center position
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
#if NDIM==3
           xc(ind,3)=(dble(iz)-0.5D0)*dx
#endif
        end do

        ! Loop over grids
        ncache=active(ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do

           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              ! Flag leaf cells
              do i=1,ngrid
                 ok(i)=son(ind_cell(i))==0
              end do

              do i=1,ngrid
                 if(ok(i))then
                    ! Get gas cell position
                    x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                    y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
#if NDIM==3
                    z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
#endif
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
#if NDIM==3
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
#else
                    dr_SN=dxx**2+dyy**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy))
#endif
                    if(((dr_SN.lt.(rSN(iSN)+0.5*dx_SN)**2) .and. (rSN(iSN) > 1.1*dx_SN)) .or. ((dr_cell < 1d-9 + dx_SN/2d0 + dx_loc/2d0) .and. (min(abs(dxx),abs(dyy),abs(dzz)) < 1d-9)))then
                       if(.not.(momentum_fb.and.(rSN(iSN) < 3d0*dx_SN)))then
                          update_boundary = .true.
                          !write(*,*)'SN blast on cpu ',myid
                          ! redistribute the mass within the SN blast uniformly and update other quantities accordingly
                          dprev=max(uold(ind_cell(i),1),smallr)
                          uold(ind_cell(i),1) = max(SNmenc(iSN)/SNvol(iSN),smallr)
!                         momprev = uold(ind_cell(i),2)**2 + uold(ind_cell(i),3)**2
!                         uold(ind_cell(i),2) = uold(ind_cell(i),2)*uold(ind_cell(i),1)/dprev
!                         uold(ind_cell(i),3) = uold(ind_cell(i),3)*uold(ind_cell(i),1)/dprev
!                         momnew = uold(ind_cell(i),2)**2 + uold(ind_cell(i),3)**2
!#if NDIM==3
!                         momprev = momprev + uold(ind_cell(i),4)**2
!                         uold(ind_cell(i),4) = uold(ind_cell(i),4)*uold(ind_cell(i),1)/dprev
!                         momnew = momnew + uold(ind_cell(i),4)**2
!#endif
!                         uold(ind_cell(i),ndim+2) = uold(ind_cell(i),ndim+2) - 0.5*momprev/dprev + 0.5*momnew/uold(ind_cell(i),1)
#if NVAR>NDIM+2+NENER
                          ! passive scalars
                          do ivar=ndim+3+nener,nvar
                             uold(ind_cell(i),ivar)=uold(ind_cell(i),ivar)*uold(ind_cell(i),1)/dprev
                          end do
#endif
!                         write(*,*)"redist on level ",ilevel,' rad ',sqrt(dr_SN),' coords ',x,' ',y,' ',z
!                         write(*,*)"redist dens:",uold(ind_cell(i),1)," redist temp:",(uold(ind_cell(i),ndim+2)*(gamma-1.0)-0.5*(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2)/uold(ind_cell(i),1))*scale_t2/uold(ind_cell(i),1)
                       endif
                       !cellweight = 1d0
!                       if ((ilevel /= level_SN(iSN)).and.(momentum_fb))then
!                          ! Calculate normalization constant for weights used to ensure correct amounts of momentum and energy injection in regions with non-uniform refinement
!                          ! Cells coarser than the central SN cell are weighted by the fraction of their volume that overlaps with the 3x3x3 finer cell SN injection region
!                          ! Cells finer than the central SN cell are weighted by how many of the fine cells that overlap with the 3x3x3 coarser cell SN injection region
!                          adjacency = 0
!                          if(abs(dxx) < 1d-9 + dx_SN/2d0)adjacency = adjacency+1
!                          if(abs(dyy) < 1d-9 + dx_SN/2d0)adjacency = adjacency+1
!                          if(abs(dzz) < 1d-9 + dx_SN/2d0)adjacency = adjacency+1
!                          if(adjacency == 0)then
!                             ! Shares a corner with central cell
!                             cellweight = 0.125d0
!                          elseif (adjacency == 1)then
!                             ! Shares an edge with central cell
!                             cellweight = 0.25d0
!                          elseif (adjacency == 2)then
!                             ! Shares a face with central cell
!                             cellweight = 0.5d0
!                          endif
!                       endif
!                       wtot(iSN) = wtot(iSN) + cellweight**(level_SN(iSN)-ilevel)
                    endif
                 endif
              end do
           end do
        end do     ! End loop over grids
     end do    ! End loop over levels
  end do  ! End loop over SNe

!  call MPI_ALLREDUCE(wtot,wtot_all,nSN,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

!  wtot=wtot_all

!  if (update_boundary)then
!    ! Update hydro quantities for split cells
!    do ilevel=nlevelmax,levelmin,-1
!       call upload_fine(ilevel)
!       do ivar=1,nvar
!          call make_virtual_fine_dp(uold(1,ivar),ilevel)
!       enddo
!    enddo
!  endif

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine subgrid_average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subgrid_Sedov_blast(xSN,mSN,rSN,vol_gas,level_SN,wtot,ncellsSN,nSN,SNfinestlevel,SNcooling,delayed)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  use sn_feedback_commons
  implicit none
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,iSN,nSN,ind,ix,iy,iz,ngrid,iskip,inds,nkin,nterm,nkin_all,nterm_all
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,ESN,vol,vol_all,dr_cell,vol_mom,vol_center
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,vol_min,dx_SN,cellweight,adjacency
  real(dp)::dr_SNs,dxxs,dyys,dzzs,cellweight_mom,cellweight_eng,xs,ys,zs,dr_cells
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_eng
  real(dp)::mom_ejecta,mom_inj,mom_term,fZ,R_cool,Tovermu,T2,nH,mu,numdens,massratio_crit,fe,nHI
  real(dp)::engfac,ektot,etherm,ektot_all,prs,fkin,R_pds,t_pds,ZonZsolar,massratio,totmom,totmom_all,einjtot,einjtot_all
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:nSN)::mSN,engdens_SN,d_gas,d_metal,vol_gas,rSN,wtot
  logical,dimension(1:nSN)::SNcooling
  real(dp),dimension(1:nSN,1:ndim)::xSN
  integer ,dimension(1:nSN)::SNfinestlevel,level_SN,ncellsSN
  logical ,dimension(1:nvector),save::ok
  logical::delayed
  integer::info

  if(nSN==0)return
  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
  vol_min=dx_min**3

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_eng = scale_d*scale_l**3*scale_v**2

  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*2d33))/scale_v**2
  do iSN=1,nSN
!     if (rSN(iSN) == 0)cycle
#ifdef DELAYED_SN
     if (delayed)then
        if (sn_isrefined(iSN)==0)cycle
     endif
#endif
     !mSN(iSN)=SN_batch_size*10d0*2d33/(scale_d*scale_l**3) ! Always assume 10 solar masses ejecta, can be changed in the future
!     d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
!     if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
     engdens_SN(iSN)=mSN(iSN)*ESN/vol_gas(iSN)
#ifdef DELAYED_SN
     if (myid==1) then
        if (SNfinestlevel(iSN) > 0) then
           write(*,*)'SN ',iSN,' will be skipped'
        else
           write(*,*)"SN at ",xSN(iSN,1), " ",xSN(iSN,2)," ",xSN(iSN,3)," ",engdens_SN(iSN)*0.67*scale_t2," ",vol_gas(iSN)
        endif
     endif
#endif
  end do

  ektot=0
  ektot_all=0
  einjtot=0
  einjtot_all=0
  totmom = 0
  totmom_all=0
  nkin=0
  nterm=0
  nkin_all = 0
  nterm_all = 0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**3
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
#if NDIM==3
        xc(ind,3)=(dble(iz)-0.5D0)*dx
#endif
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
#if NDIM==3
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
#endif
                 do iSN=1,nSN
#ifdef DELAYED_SN
                    if (delayed)then
                       if (sn_isrefined(iSN)==0)cycle
                    endif
                    if (SNfinestlevel(iSN) == 0)then
#endif
                    if(rSN(iSN)<0)cycle
                       dx_SN = scale*0.5D0**level_SN(iSN)
                       vol_center = dx_SN**3
                       vol_mom = vol_gas(iSN) - vol_center
                       ! Check if the cell lies within the SN radius
                       dxx=x-xSN(iSN,1)
                       dyy=y-xSN(iSN,2)
#if NDIM==3
                       dzz=z-xSN(iSN,3)
                       dr_SN=dxx**2+dyy**2+dzz**2
                       dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
#else
                       dr_SN=dxx**2+dyy**2
                       dr_cell=MAX(ABS(dxx),ABS(dyy))
#endif
                       if (rSN(iSN) == 0)then
                          if(dr_SN < 1d-10)uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2) + engdens_SN(iSN)
                       elseif(((dr_SN.lt.(rSN(iSN)+0.5*dx_SN)**2) .and. (rSN(iSN) > 1.1*dx_SN)) .or. ((dr_cell < 1d-9 + dx_SN/2d0 + dx_loc/2d0) .and. (min(abs(dxx),abs(dyy),abs(dzz)) < 1d-9)))then
                          if (momentum_fb)then
                              ! Kinetic feedback: Inject some fraction of SN energy as kinetic energy to alleviate overcooling
                              ! This largely follows either Gentry, Madau & Krumholz (2020) or Simpson et al. (2015)
                              ! but without star particles and no mass or metals is injected
                              mom_ejecta = (vol_mom/vol_gas(iSN))*sqrt(2d0*mSN(iSN)*1d51/scale_eng)
                              engfac = 1d0
                              cellweight = 1d0
                              if (dr_SN > 1d-10)then
                                 dr_SN = sqrt(dr_SN)
                                   adjacency = 0
                                   if(abs(dxx) < 1d-9)adjacency = adjacency+1
                                   if(abs(dyy) < 1d-9)adjacency = adjacency+1
                                   if(abs(dzz) < 1d-9)adjacency = adjacency+1
                                   !write(*,*)"dxx,dyy,dzz,adj",dxx,dyy,dzz,adjacency
                                   cellweight_mom = (adjacency/2d0)*1.5d0!0.125d0

!                                 if (ilevel /= level_SN(iSN))then
!                                   adjacency = 0
!                                   if(abs(dxx) < 1d-9 + dx_SN/2d0)adjacency = adjacency+1
!                                   if(abs(dyy) < 1d-9 + dx_SN/2d0)adjacency = adjacency+1
!                                   if(abs(dzz) < 1d-9 + dx_SN/2d0)adjacency = adjacency+1
!                                   write(*,*)"dxx,dyy,dzz,adj",dxx,dyy,dzz,adjacency
!                                   if(adjacency == 0)then
!                                      ! Shares a corner with central cell
!                                      cellweight = 0.125d0
!                                   elseif (adjacency == 1)then
!                                      ! Shares an edge with central cell
!                                      cellweight = 0.25d0
!                                   elseif (adjacency == 2)then
!                                      ! Shares a face with central cell
!                                      cellweight = 0.5d0
!                                   endif
!!                                   cellweight = 0.125d0
!                                 endif
!                                 ! Momentum injection region excludes central cell and so is weighted differently from energy and mass
!                                 cellweight_mom = (cellweight**(level_SN(iSN)-ilevel))!*(ncellsSN(iSN)-1)/(wtot(iSN)-1d0)
!!!                                 cellweight_mom = 1d0
!                                 if (ilevel < level_SN(iSN))cellweight_mom = 0.125d0
!                                 if (ilevel > level_SN(iSN))cellweight_mom = 8d0
                                 cellweight_eng = (adjacency/2d0)*(19d0/13d0)!cellweight_mom!(cellweight**(level_SN(iSN)-ilevel))!*ncellsSN(iSN)/wtot(iSN)
                                 massratio = sqrt(max(uold(ind_cell(i),1),smallr)*vol_center/(cellweight_eng*mSN(iSN)/ncellsSN(iSN))) !vol_gas(iSN)/(mSN(iSN)))
                                 prs = (uold(ind_cell(i),ndim+2) - 0.5d0*(uold(ind_cell(i),2)**2 + uold(ind_cell(i),3)**2 + uold(ind_cell(i),4)**2)/max(uold(ind_cell(i),1),smallr))*(gamma-1.0)
                                 Tovermu = prs/max(uold(ind_cell(i),1),smallr)*scale_T2
                                 nH = max(uold(ind_cell(i),1),smallr)*scale_nH
                                 call GetMuAndTemperature(Tovermu,nH,mu,T2,nHI)
                                 numdens = max(uold(ind_cell(i),1),smallr)/mu
                                 ZonZsolar = uold(ind_cell(i),imetal)/uold(ind_cell(i),1)/0.02
                                 if (simpson_fb)then
                                    ! Use scheme of Simpson et al. (2015) to calculate fraction of kinetic energy
                                    if (ZonZsolar > 0.01)then
                                       t_pds = 26.5d0 * numdens**(-4d0/7d0)*ZonZsolar**(-5d0/14d0) ! in kyr
                                       R_pds = 18.5d0 * numdens**(-3d0/7d0)*ZonZsolar**(-1d0/7d0)  ! in pc
                                    else
                                       t_pds = 3.06d2 * numdens**(-3d0/4d0) ! in kyr
                                       R_pds = 49.3d0 * numdens**(-1d0/2d0) ! in pc
                                    endif
                                    if (R_pds > 4.5d0*1d3*dx_loc)then
                                       fkin = 0.0 ! Reduce to thermal dump
                                    else
                                       fkin = 3.97d-6 * max(uold(ind_cell(i),1),smallr)*R_pds**7*t_pds**(-2)*(1d3*dx_loc)**(-2)
                                       if (fkin > 1d0)fkin = 1d0
                                    endif
                                    mom_inj = cellweight_mom*sqrt(fkin)*mom_ejecta*massratio/vol_mom
                                    write(*,*)"fkin: ",fkin," t_pds (kyr): ",t_pds," R_pds (kpc): ",1d-3*R_pds
                                 else
                                    ! Use scheme of Gentry, Madau & Krumholz (2020) to inject either terminal momentum or 100% kinetic energy
                                    fZ = 2d0
                                    if (ZonZsolar > 0.01)fZ=ZonZsolar**(-0.14)
                                    if (cioffi_mom)then
                                       mom_term = mom_fac * 9.6d43 * SN_batch_size**(13d0/14d0)*numdens**(-1d0/7d0)*fZ**(3d0/2d0)/(scale_d*scale_l**3*scale_v) !Terminal momentum from Cioffi+ 1988
                                    else
                                       mom_term = mom_fac * 6d43 * SN_batch_size**(16d0/17d0)*numdens**(-2d0/17d0)*fZ/(scale_d*scale_l**3*scale_v) !Terminal momentum from Rosdahl+ 2017
                                    endif
                                    !if (sqrt(1d0 + max(uold(ind_cell(i),1),smallr)*vol_gas(iSN)/mSN(iSN)) < mom_term/mom_ejecta)then
                                    !   mom_inj = mom_ejecta*massratio/vol_gas(iSN)
                                    !else
                                    !   mom_inj = mom_term/vol_gas(iSN)
                                    !endif
                                    if (massratio > mom_term/mom_ejecta)then
                                       nterm=nterm+1
                                    else
                                       nkin=nkin+1
                                    endif
                                    if (sn_smooth_transition)then
                                       massratio_crit = 900d0/(mSN(iSN)*(scale_d*scale_l**3/2d33)*0.667d0)*SN_batch_size**(-2d0/17d0)*numdens**(-4d0/17d0)*ZonZsolar**(-0.28d0)
                                       fe = 1d0 - 0.333d0*(massratio**2d0 - 1d0)/(massratio_crit - 1d0)
                                       if (massratio**2d0 < massratio_crit)then
                                          mom_inj = cellweight_mom*mom_ejecta*massratio*sqrt(fe)/vol_mom
                                       else
                                          mom_inj = cellweight_mom*mom_term/vol_mom
                                       endif
                                    else
                                       mom_inj = cellweight_mom*mom_ejecta*min(massratio, mom_term/mom_ejecta)/vol_mom
                                    endif
                                 endif
                                 if ((ilevel < level_SN(iSN)) .and. (adjacency > 0))then
                                    ! Sample subcells of coarse cell overlapping 2 or 4 SN injection region cells on finer level
                                    do inds=1,twotondim
                                       iz=(inds-1)/4
                                       iy=(inds-1-4*iz)/2
                                       ix=(inds-1-2*iy-4*iz)
                                       xs=x+(dble(ix)-0.5D0)*dx_SN
                                       ys=y+(dble(iy)-0.5D0)*dx_SN
                                       dxxs=xs-xSN(iSN,1)
                                       dyys=ys-xSN(iSN,2)
#if NDIM==3
                                       zs=z+(dble(iz)-0.5D0)*dx_SN
                                       dzzs=zs-xSN(iSN,3)
                                       dr_SNs=dxxs**2+dyys**2+dzzs**2
                                       dr_cells=MAX(ABS(dxxs),ABS(dyys),ABS(dzzs))
#else
                                       dr_SNs=dxxs**2+dyys**2
                                       dr_cells=MAX(ABS(dxxs),ABS(dyys))
#endif
                                       dr_SNs = sqrt(dr_SNs)
                                       if(dr_cells < 1d-9 + dx_SN)then
                                          uold(ind_cell(i),2)=uold(ind_cell(i),2) + mom_inj*dxxs/dr_SNs
                                          uold(ind_cell(i),3)=uold(ind_cell(i),3) + mom_inj*dyys/dr_SNs
#if NDIM==3
                                          uold(ind_cell(i),4)=uold(ind_cell(i),4) + mom_inj*dzzs/dr_SNs
#endif
                                          uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2) + cellweight_eng*engdens_SN(iSN)
                                          write(*,*)"injection in subcell: ",xs,ys,zs,dxxs,dyys,dzzs,ix,iy,iz,8.0*mom_inj*dyys/dr_SNs,dr_SNs
                                       endif
                                    enddo
                                 else
                                    uold(ind_cell(i),2)=uold(ind_cell(i),2) + mom_inj*dxx/dr_SN
                                    uold(ind_cell(i),3)=uold(ind_cell(i),3) + mom_inj*dyy/dr_SN
#if NDIM==3
                                    uold(ind_cell(i),4)=uold(ind_cell(i),4) + mom_inj*dzz/dr_SN
#endif
                                    R_cool = 0.0284*SN_batch_size**(2d0/7d0)*numdens**(-3d0/7d0)*fZ
                                    if ((dr_SN > R_cool).and.(ilevel >= level_SN(ilevel)).and.Rcool_correction)then
                                       engfac = (dr_SN/R_cool)**(-6.5d0)
                                       etherm = (engdens_SN(iSN) - 0.5d0*((mom_inj*dxx/dr_SN)**2 + (mom_inj*dyy/dr_SN)**2 + (mom_inj*dzz/dr_SN)**2)/uold(ind_cell(i),1))*engfac
                                       engdens_SN(iSN) = etherm + 0.5d0*((mom_inj*dxx/dr_SN)**2 + (mom_inj*dyy/dr_SN)**2 + (mom_inj*dzz/dr_SN)**2)/uold(ind_cell(i),1)
                                    endif
                                    uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2) + cellweight_eng*engdens_SN(iSN)
                                 endif
                                 if(index(output_dir,"sntest") > 0)then
                                    write(*,*)"Tovermu, T, mu, numdens, vol_gas, Rcool, engfac, mom_inj, mom_term, e_inj: ",Tovermu, T2, mu, numdens, vol_gas(iSN), R_cool, engfac, mom_inj*vol_gas(iSN), mom_term, engdens_SN(iSN)
                                    write(*,*)"dx_loc", dx_loc,"w_mom",cellweight_mom,"w_eng",cellweight_eng,"wtot",wtot(iSN),"x ",x," y ",y," z ",z," dxx",dxx," dyy",dyy," dzz",dzz," dr_SN",dr_SN," momx",mom_inj*dxx/dr_SN,"momy ",mom_inj*dyy/dr_SN," momz ",mom_inj*dzz/dr_SN
                                 endif
                              else
                                 cellweight_eng=(19d0/13d0)!ncellsSN(iSN)/wtot(iSN)
                                 uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2) + cellweight_eng*engdens_SN(iSN)
                              endif
!                              uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2) + cellweight_eng*engdens_SN(iSN)
                          else
                              ! Thermal dump
                              ! Update the total energy of the gas
                              uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+engdens_SN(iSN)
                              if (.not.SNcooling(iSN))uold(ind_cell(i),idelay) = uold(ind_cell(i),idelay) + max(uold(ind_cell(i),1),smallr) ! + d_gas(iSN)
                          endif
                          etherm = uold(ind_cell(i),ndim+2) - 0.5d0*(uold(ind_cell(i),2)**2 + uold(ind_cell(i),3)**2 + uold(ind_cell(i),4)**2)/uold(ind_cell(i),1)
                          ektot = ektot + (uold(ind_cell(i),ndim+2) - etherm)*vol_loc
                          einjtot = einjtot + cellweight_eng*engdens_SN(iSN)*vol_loc
                          totmom=totmom+(abs(uold(ind_cell(i),2)) + abs(uold(ind_cell(i),3)) + abs(uold(ind_cell(i),4)))*vol_center
                       endif
#ifdef DELAYED_SN
                   endif
#endif
                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels
  
#ifndef WITHOUTMPI
  call MPI_REDUCE(ektot,ektot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
  call MPI_REDUCE(einjtot,einjtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
  call MPI_REDUCE(totmom,totmom_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
  call MPI_REDUCE(nkin,nkin_all,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,info)
  call MPI_REDUCE(nterm,nterm_all,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,info)
  if (myid==1)write(*,*)"Ekin: ",ektot_all*scale_eng," Einj: ", einjtot_all*scale_eng, " mom:",totmom_all*(scale_d*scale_l**3*scale_v/(2e33*1e5))
  if (myid==1)write(*,*)"N_kin: ",nkin_all," N_term:",nterm_all
!  if (myid==1)write(*,*)"Ekin: ",ektot_all*scale_eng
#endif

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine subgrid_Sedov_blast


! Routines to output and read seeds used for the supernova related random number generation
! Ensures that the SN feedback continues properly when restarting a simulation
subroutine output_seeds(filename)
  use amr_commons
  use hydro_commons
  use mpi_mod
  use sn_feedback_commons
  implicit none

  integer::info
  integer ,dimension(1:IRandNumSize,1:ncpu)::allseeds
  character(LEN=256)::filename, fileloc
  integer::ilun,icpu

  if(verbose)write(*,*)'Entering output_seeds'
  if(myid==1)then
     ! Open file
     fileloc=TRIM(filename)
     open(newunit=ilun,file=fileloc,form='formatted')
     write(ilun,*)"ProcID  Seeds"
  end if
  
  call MPI_GATHER(localseedsn,IRandNumSize,MPI_INTEGER,allseeds,IRandNumSize,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  if(myid==1)then
      do icpu=1,ncpu
         write(ilun,'(I6,4I24)') icpu, allseeds(1,icpu), allseeds(2,icpu), allseeds(3,icpu), allseeds(4,icpu) ! Assumes IRandNumSize=4
      enddo
      close(ilun)
  endif
  if(verbose)write(*,*)'Exiting output_seeds'
end subroutine output_seeds


subroutine init_subgrid_feedback
  use amr_commons
  use hydro_commons
  use mpi_mod
  use sn_feedback_commons
  implicit none
  integer::info
  integer ,dimension(1:IRandNumSize,1:ncpu)::allseeds
  character(LEN=256)::fileloc,dummy
  character(LEN=5)::nchar
  integer::ilun,icpu,ccpu
  logical::file_exist

  if(verbose)write(*,*)'Entering init_subgrid_feedback'
  if(nrestart>0)then
      call title(nrestart,nchar)
      fileloc=trim(output_dir)//'output_'//TRIM(nchar)//'/subgrid_sn_seeds'//TRIM(nchar)//'.txt'
      inquire(file=fileloc,exist=file_exist)
      if(file_exist) then
         if(myid==1)then
            open(newunit=ilun,file=fileloc)
            read(ilun,*)dummy
            do icpu=1,ncpu
               read(ilun,*) ccpu, allseeds(:,icpu)
            enddo
            close(ilun)
         endif
         call MPI_SCATTER(allseeds,IRandNumSize,MPI_INTEGER,localseedsn,IRandNumSize,MPI_INTEGER,0,MPI_COMM_WORLD,info)
      else
         iseed=seed_init
      endif
  else
     iseed=seed_init
  endif

  if (myid==200) then
     write(*,*)'seeds(200): ',localseedsn
  endif
  if(verbose)write(*,*)'Exiting init_subgrid_feedback'
end subroutine init_subgrid_feedback


subroutine GetMuAndTemperature(Tovermu,nH,mu,T2,nHI)
  use amr_commons, ONLY: t
  use amr_parameters, ONLY: dp, aexp, evolve_uvb
  use cooling_module, ONLY: set_rates, cmp_chem_eq, get_uvb_expfac
  implicit none
  real(dp)::Tovermu,T2,nH,mu,nHI,z,cura
  real(dp)::mu_old,err_mu,mu_left,mu_right,n_TOT
  real(dp),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(dp),dimension(1:6) :: n_spec
  integer::niter

  if(evolve_uvb)then
     cura = get_uvb_expfac(t)
  else
     cura = aexp
  endif
  call set_rates(t_rad_spec,h_rad_spec,cura)
  z = 1.d0/cura-1.D0

  ! Iteration to find mu
  err_mu=1.
  mu_left=0.5
  mu_right=1.3
  niter=0
  do while (err_mu > 1.d-4 .and. niter <= 50)
     mu_old=0.5*(mu_left+mu_right)
     T2 = Tovermu*mu_old
     call cmp_chem_eq(T2,nH,t_rad_spec,n_spec,n_TOT,mu,z)
     nHI = n_spec(2)
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
end subroutine GetMuAndTemperature

!###########################################################
!###########################################################
!###########################################################
!###########################################################

