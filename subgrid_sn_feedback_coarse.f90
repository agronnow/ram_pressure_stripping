!------------------------------------------------------------------------------------------------------------
! Supernova feedback without star particles by Asger Gronnow
! NOTE: This version should only be called on levelmin, use only for low SFR galaxies
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
subroutine subgrid_sn_feedback_coarse(icount)
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

!  if (ilevel/=levelmin)return

  if(numbtot(1,levelmin)==0) return
  if(.not. hydro)return
  !if(ndim.ne.3)return
  if(static)return

  if(verbose)write(*,*)' Entering SN feedback'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
!  dx=0.5D0**ilevel
!  nx_loc=(icoarse_max-icoarse_min+1)
!  skip_loc=(/0.0d0,0.0d0,0.0d0/)
!  if(ndim>0)skip_loc(1)=dble(icoarse_min)
!  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
!  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
!  scale=boxlen/dble(nx_loc)
!  dx_loc=dx*scale
!  vol_loc=dx_loc**3
!  dx_min=(0.5D0**nlevelmax)*scale
!  vol_min=dx_min**3

  ESN=1d51/(10.*2d33)/scale_v**2

  scale_m = scale_d*scale_l**3

  ! Cells center position relative to grid center position
!  do ind=1,twotondim
!     iz=(ind-1)/4
!     iy=(ind-1-4*iz)/2
!     ix=(ind-1-2*iy-4*iz)
!     xc(ind,1)=(dble(ix)-0.5D0)*dx
!     xc(ind,2)=(dble(iy)-0.5D0)*dx
!     xc(ind,3)=(dble(iz)-0.5D0)*dx
!  end do

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

  doSNIa = .true.!.false.
!#if defined(SNIA_FEEDBACK) && !defined(SN_INJECT)
!  if (ilevel == nlevelmax) then
!     doSNIa = .true.
! Only generate SNIa at the last subcycle of the finest level
!     if (icount == 2) then
!        if(nsubcycle(ilevel)==2)doSNIa = .true.
!     else
!        doSNIa = .true.
!     endif
!  endif
!#endif

  nSNIa = 0
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
  
  if ((myid==1) .and. ((pot_grow_rate > 0.0) .or. firstcall) .and. .not.(pot_rad_max)) then
    ! Calculate SNIa radius Comulative Distribution Function for the current potential truncation radius
    potrad = r_cut*(1d0+pot_grow_rate*(t-t_pot_grow_start))
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
#ifdef DEBUG_SNIA
    if ((t==0.0) .or. pot_rad_max)then
       fileloc=trim(output_dir)//sniadist_fname
       ilun=130
       open(ilun, file=fileloc, form='formatted')
       write(ilun,*)"rad                      PDF                   CDF"
       do i=1,NPDFBINS
          write(ilun,'(3E26.16)') (i-1)*binwidth, PDF_SNIa(i), CDF_SNIa(i)
       enddo
       close(ilun)
       write(*,*)"Wrote ",sniadist_fname
    endif
#endif
  endif

  if (.not.(firstcall)) then
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
    if(verbose.and.(myid==1))write(*,*)"Calculating SNIa"
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
      PoissMeanIa=rho_SNIa*scale_t/(3600*24*365.25)*dtnew(levelmin)/SN_batch_size ! Get expected number of SNIa formed taking units into account
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
      allocate(min_r2_all(1:2,1:nSNIa))
      call MPI_ALLREDUCE(min_r2,min_r2_all,nSNIa,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,info)
    endif
  endif
    
#endif

  !------------------------------------------------
  ! Compute number of SNe in each cell
  !------------------------------------------------

#ifdef SNIA_FEEDBACK
  sfr = 0.0
#endif

  pdevsn2=0
  nSN_loc = 0

  ! Mesh spacing in that level                                                                                                                                                                                                                                                            
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

do ilevel=levelmin,nlevelmax
  if(verbose.and.(myid==1))write(*,*)"Calculating SN on level",ilevel
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  dx_loc=dx*scale
  vol_loc=dx_loc**3
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**3

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
! get values of uold for density and velocities in virtual boundaries
!#ifndef WITHOUTMPI
!  do ivar=1,nvar
!     call make_virtual_fine_dp(uold(1,ivar),ilevel)
!  end do
!#endif

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
                ! Poisson mean
                PoissMean=rho_SN*rho_sfr*scale_t/(3600*24*365.25)*vol_loc*dtnew(levelmin)/SN_batch_size ! Get expected number of SNe formed taking units into account
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
                   nSN_loc=nSN_loc+1
                   xSN_loc(nSN_loc,1)=x
                   xSN_loc(nSN_loc,2)=y
#if NDIM==3
                   xSN_loc(nSN_loc,3)=z
#endif
                   levelSN_loc(nSN_loc) = ilevel
                   mSN_loc(nSN_loc)=SN_batch_size*10.0*2d33/(scale_d*scale_l**3) !Always assume 10 solar mass ejection
                   !if (outputSN) then
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
                   !endif
                !   do iSN=1,nSN
                !      uold(ind_cell(i),ndim+2) = uold(ind_cell(i),ndim+2) + 10.0*2d33/(scale_d*scale_l**3)*ESN/vol_loc
                !      write(*,*) "SN explosion!"
                !   end do
                end do
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

     do iSN=1,nSNIa
        if ((int(min_r2_all(2,iSN)) == myid) .and. (levelSNIa(iSN) == ilevel)) then ! The cell containing the SN center is on this processor and level
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
           write(ilun,'(4E26.16,E26.16,2I5)') t, xSNIa(iSN,1), xSNIa(iSN,2), z, PoissMeanIa, levelSNIa(iSN),myid!, oldseed(1), oldseed(2), oldseed(3), oldseed(4)
           close(ilun)
           !endif
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
  sfr_tot(ilevel-levelmin+1) = sfr_tot(ilevel-levelmin+1) + sfr_tot_level
enddo !loop over levels

  if(verbose.and.(myid==1))write(*,*)"Done SN loop over levels"
  if (calc_sfr) then
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
  if(verbose.and.(myid==1))write(*,*)"Tot SN ",nSN_tot
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

  if(verbose.and.(myid==1))write(*,*)"Applying SN"
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
  if (nSNIa>0)deallocate(xpdf,xSNIa,min_r2,min_r2_all,levelSNIa)
#endif

  ! Update hydro quantities for split cells
!  do ilevel=nlevelmax,levelmin,-1
!     call upload_fine(ilevel)
!     do ivar=1,nvar
!        call make_virtual_fine_dp(uold(1,ivar),ilevel)
!     enddo
!  enddo

  if(verbose)write(*,*)'Exiting subgrid_sn_feedback_coarse'

end subroutine subgrid_sn_feedback_coarse
!#endif

