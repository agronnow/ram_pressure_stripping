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
  integer,parameter::RADCELL_MAX=10


#ifndef WITHOUTMPI
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all
  real(dp),dimension(:,:),allocatable::xSN_all
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
  integer,dimension(:),allocatable::indSN
  real(dp),dimension(:),allocatable::mSN,mSN_loc,ekBlast,rSN
  real(dp),dimension(:,:),allocatable::xSN,xSN_loc,dq,vol_gas

  logical,save::firstcall = .true.
  logical,save::calc_sfr = .false.
  real(dp),save::t_sfrlog = 0.0
  real(dp)::cursfr = 0.0
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
  real(dp)::unif_rand,r2,rho_SNIa,area,DTD_A,DTD_s,currad,Phi0,PhiR,mu_cloud,rho0dm,r,Pinf
  real(dp)::PoissMeanIa,c_s2,nHc,ctime,csfh,diff,mindiff,signx,dt,sfr,sfr_tot_level
  integer::nSNIa,nt,imin,stat
  real(dp),dimension(:,:),allocatable::xpdf,xSNIa,min_r2,min_r2_all
  logical::doSNIa
#endif

  logical,save::firstsn = .true.

  ! TODO: when f2008 is obligatory - remove this and replace erfc_pre_f08 below by
  ! the f2008 intrinsic erfc() function:
  real(dp) erfc_pre_f08

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

  ! Star formation time scale from Gyr to code units
  ! SFR apply here for long lived stars only
!  t0=t_star*(1d9*365.*24.*3600.)/scale_t
!  trel=sf_trelax*1d6*(365.*24.*3600.)/scale_t

  ! ISM density threshold from H/cc to code units
!  nISM = n_star
!  if(cosmo)then
!     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
!     nISM = MAX(nCOM,nISM)
!  endif
!  d0   = nISM/scale_nH
!  d00  = n_star/scale_nH

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

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

!  !------------------------------------------------
!  ! Convert hydro variables to primitive variables
!  !------------------------------------------------
!  ncache=active(ilevel)%ngrid
!  do igrid=1,ncache,nvector
!     ngrid=MIN(nvector,ncache-igrid+1)
!     do i=1,ngrid
!        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
!     end do
!     do ind=1,twotondim
!        iskip=ncoarse+(ind-1)*ngridmax
!        do i=1,ngrid
!           ind_cell(i)=iskip+ind_grid(i)
!        end do
!        do i=1,ngrid
!           d=uold(ind_cell(i),1)
!           u=uold(ind_cell(i),2)/d
!           v=uold(ind_cell(i),3)/d
!#if NDIM==3
!!           w=uold(ind_cell(i),4)/d
!#endif
!           e=uold(ind_cell(i),ndim+2)
!#ifdef SOLVERmhd
!           bx1=uold(ind_cell(i),6)
!           by1=uold(ind_cell(i),7)
!#if NDIM==3
!           bz1=uold(ind_cell(i),8)
!#endif
!           bx2=uold(ind_cell(i),nvar+1)
!           by2=uold(ind_cell(i),nvar+2)
!#if NDIM==3
!           bz2=uold(ind_cell(i),nvar+3)
!           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
!#else
!           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2)
!#endif
!#endif
!           e=e-0.5d0*d*(u**2+v**2+w**2)
!#if NENER>0
!           do irad=0,nener-1
!              e=e-uold(ind_cell(i),inener+irad)
!           end do
!#endif
!           uold(ind_cell(i),1)=d
!           uold(ind_cell(i),2)=u
!           uold(ind_cell(i),3)=v
!#if NDIM==3
!           uold(ind_cell(i),4)=w
!#endif
!           prs = e/d
!           !uold(ind_cell(i),5)=e/d
!        end do
!        do ivar=imetal,nvar
!           do i=1,ngrid
!              d=uold(ind_cell(i),1)
!              w=uold(ind_cell(i),ivar)/d
!              uold(ind_cell(i),ivar)=w
!           end do
!        end do
!     end do
!  end do

! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  allocate(xSN_loc(1:100,1:ndim),mSN_loc(1:100))

#ifdef SNIA_FEEDBACK
#ifdef SNIA_HERINGER19
  DTD_A = 10.0**(-12.15) ! +0.10 dex -0.13 dex, [M_sun^-1 yr^-1], from Heringer et al. (2019)
  DTD_s = -1.34          ! +0.19 -0.17, from Heringer et al. (2019)
#else
  DTD_A = 2.11d-13       ! [M_sun^-1 yr^-1], from Maoz et al. (2012)
  DTD_s = -1.12          ! From Maoz et al. (2012)
#endif

  doSNIa = .true.
!  doSNIa = .false.
!  if (ilevel == nlevelmax) then ! Only generate SNIa at the last subcycle of the finest level
!     if (icount == 2) then
!        if(nsubcycle(ilevel)==2)doSNIa = .true.
!     else
!        doSNIa = .true.
!     endif
!  endif

  nSNIa = 0
  if ((ilevel == levelmin) .and. (calc_sfr)) then
    calc_sfr = .false.
    if (sfhist_update) then
       nhist = nhist + 1
       t_sfhist(nhist) = t*scale_t/3.154e16 + tinit_sim
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
     fileloc=trim(sfhistfile)
     ilun=150
     inquire(file=fileloc,exist=file_exist)
     nhist = 0
     if(file_exist) then
        open(ilun, file=fileloc)
        read(ilun,*)dummyline
        do
           read(ilun,*, iostat=stat)ctime,csfh
           if ((stat /= 0) .or. (ctime > tinit_sim))exit
           t_sfhist(nhist+1) = ctime
           if (t_sfhist(nhist+1) < 0.1)t_sfhist(nhist+1) = 0.1
           sfhist(nhist+1) = csfh
           nhist=nhist+1
           if (myid==1)write(*,*)nhist,t_sfhist(nhist),sfhist(nhist)
        end do
!        t_sfhfile = t_sfhist(nhist)
        close(ilun)
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
           if ((stat /= 0) .or. (ctime > tinit_sim + t*scale_t/3.154e16))exit
           t_sfhist(nhist+1) = ctime
           if (t_sfhist(nhist+1) < 0.1)t_sfhist(nhist+1) = 0.1
           sfhist(nhist+1) = csfh
           nhist=nhist+1
           write(*,*)myid,nhist,t_sfhist(nhist),sfhist(nhist)
        end do
!        t_sfhfile = t_sfhist(nhist)
        close(ilun)
     endif

     if (myid==1) then
       ! Calculate SNIa radius Comulative Distribution Function
#ifndef SNIA_PLUMMER
       nHc = n0g*scale_nH
       call GetMuFromTemperature(T_cloud,nHc,mu_cloud)
       rho0dm = gravity_params(1)
       Phi0 = -2.0*twopi*rho0dm*R_s**2
       c_s2 = kb*T_cloud/(mu_cloud*mh)/scale_v**2 !square of isothermal sound speed in cloud centre
       Pinf= dexp(-(Phi0*R_s*dlog(1+Rad_cloud/R_s)/Rad_cloud - Phi0)/c_s2)
#endif
       area = 0.0
       binwidth = Rad_cloud/(NPDFBINS-1.0)
       
       do i=1,NPDFBINS
          currad = (i-1)*binwidth
#ifdef SNIA_PLUMMER
          PDF_SNIa(i) = (3.0/(2.0*twopi*r_plummer**3))*(1.0+currad/r_plummer**2)**(-2.5)
#else
          if (currad > 0.0) then
            PhiR = Phi0*R_s*dlog(1+currad/R_s)/currad
          else
            PhiR = Phi0
          endif
          PDF_SNIa(i) = dexp(-(PhiR-Phi0)/c_s2)-Pinf
#endif
          area = area + PDF_SNIa(i)*binwidth
       enddo
       do i=1,NPDFBINS
          ! Normalize PDF
          PDF_SNIa(i) = PDF_SNIa(i)/area
       enddo
       do i=1,NPDFBINS
          CDF_SNIa(i) = sum(PDF_SNIa(1:i))*binwidth
       enddo
       xCloud(1) = x1_c
       xCloud(2) = x2_c
#if NDIM==3
       xCloud(3) = x3_c
#endif

#ifdef DEBUG_SNIA
       fileloc=trim(output_dir)//'snIa_pdf.dat'
       ilun=130
       open(ilun, file=fileloc, form='formatted')
       write(ilun,*)"rad                      PDF                   CDF"
       do i=1,NPDFBINS
          write(ilun,'(3E26.16)') (i-1)*binwidth, PDF_SNIa(i), CDF_SNIa(i)
       enddo
       close(ilun)
#endif
     endif
  endif

  if ((ilevel == levelmin) .and. .not.(firstcall)) then
    if (t*scale_t/3.154e16 + tinit_sim - t_sfhist(nhist) > dt_sfhist) then
       if (myid==1)write(*,*)'Update sfh',t*scale_t/3.154e16,tinit_sim,t_sfhist(nhist),dt_sfhist
       sfhist_update = .true.
       calc_sfr = .true.
    else if (t*scale_t/3.154e16 + tinit_sim - t_sfrlog > dt_sfrlog) then
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
            dt = t*scale_t/3.154e16 + tinit_sim - t_sfhist(nhist)
         endif
         rho_SNIa = rho_SNIa + sfhist(nhist-nt+1)*DTD_A*(t_sfhist(nt))**DTD_s*dt*1d9 ! SNIa per year
         !write(*,*)nt,nhist-nt+1,rho_SNIa,t_sfhist(nt),sfhist(nhist-nt+1),dt,DTD_A,DTD_s
      end do
      PoissMeanIa=rho_SNIa*scale_t/(3600*24*365.25)*dtnew(ilevel) ! Get expected number of SNIa formed taking units into account
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
      write(ilun,'(3E26.16)') t, rho_SNIa, PoissMeanIa
      close(ilun)
#endif
      ! Compute Poisson realisation
      call poissdev(localseed,PoissMeanIa,nSNIa)
      if (nSNIa > 0) then
         allocate(xpdf(1:nSNIa,1:ndim),xSNIa(1:nSNIa,1:ndim))
         ! Generate random positions for SNIa assuming that the stellar distribution follows the initial distribution of star forming gas
         do iSN=1,nSNIa
            do idim=1,ndim
               ! Generate random real on [0,1], transform to [-1,1], store the sign and transform back to [0,1] to allow coordinates on either side of the cloud
               call ranf(localseed, unif_rand)
               r = 2.0*unif_rand - 1.0
               signx = sign(1d0,r)
               r = abs(r)
               mindiff = 1.e22
               imin = 1
               do i=1,NPDFBINS
                  diff = abs(CDF_SNIa(i) - r)
                  if (diff < mindiff) then
                     mindiff = diff
                     imin = i
                  endif
               enddo
               xpdf(iSN,idim) = signx*binwidth*(imin-1) + xCloud(idim)*boxlen ! Random coordinate along each axis
write(*,*)'nSNIa',nSNIa,'SNIa',iSN,'unif_rand',unif_rand,'signx',signx,'r',r,'x(',idim,')=',xpdf(iSN,idim)
            enddo
         enddo
      endif
    endif

    call MPI_BCAST(nSNIa,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    if (nSNIa > 0) then
      if (myid/=1)allocate(xpdf(1:nSNIa,1:ndim),xSNIa(1:nSNIa,1:ndim))
      allocate(min_r2(1:2,1:nSNIa))
      call MPI_BARRIER(MPI_COMM_WORLD,info)
      call MPI_BCAST(xpdf,nSNIa*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
      call MPI_BCAST(xSNIa,nSNIa*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
      call MPI_BCAST(PoissMeanIa,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)

      do iSN=1,nSNIa
        min_r2(1,iSN) = 1.e22 ! Huge value, must always be greater than min distance to SN
        min_r2(2,iSN) = myid  ! Tag with rank on each cpu for later use with MPI_MINLOC
      enddo
      write(*,*) 'SNIa explosion'
    endif
  endif
    
#endif

  !------------------------------------------------
  ! Compute number of SNe in each cell
  !------------------------------------------------
!  ntot=0
!  ndebris_tot=0
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
           d=uold(ind_cell(i),1)
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
           nH=max(uold(ind_cell(i),1),smallr)*scale_nH
!          T_poly=T2_star*(nH/nISM)**(g_star-1.0)
!          T2=T2-T_poly
           if(T2>4e4) then
              ok(i)=.false.
              Tfail=Tfail+1
           else
              Tpass=Tpass+1
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
        ! Calculate number of new stars in each cell using Poisson statistics
        !do i=1,ngrid
        !   nstar(i)=0
           ! Poisson mean
        !   rho_sfr = 10.0**(0.9+1.91*dlog10(uold(ind_cell(i),1)*0.02439303604))
        !   PoissMean=rho_sfr*scale_t/(3600*24*365.25)*vol_loc*dtnew(ilevel)*2d33/scale_m/mstar
        !   tot_sf = tot_sf + PoissMean
        !enddo

#ifdef SNIA_FEEDBACK
        if (doSNIa) then
          ! Find location of SNIa, i.e. the cell closest to the randomly generated position
          do iSN=1,nSNIa
            do i=1,ngrid
              if (son(ind_cell(i))==0) then
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
                 endif
              endif
            enddo
          enddo
        endif
#endif

        !SNII
        do i=1,ngrid
           if(ok(i))then !Compute number of SNII if temperature is low enough for star formation
              d=uold(ind_cell(i),imetal+1) !SF gas density (gas initially within r_SFR)
              !mcell=d*vol_loc
              ! Free fall time of an homogeneous sphere
              !tstar= .5427*sqrt(1.0/(factG*max(d,smallr)))
              !if(.not.sf_virial) sfr_ff(i) = eps_star
              ! Gas mass to be converted into stars
              !mgas=dtnew(ilevel)*(sfr_ff(i)/tstar)*mcell
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
                PoissMean=rho_SN*rho_sfr*scale_t/(3600*24*365.25)*vol_loc*dtnew(ilevel) ! Get expected number of SNe formed taking units into account
                !PoissMean=mgas/mstar
                !write(*,*) "PoissMean: ", PoissMean, " rho_sfr: ", rho_sfr, " dt: ", dtnew(ilevel), " vol: ", vol_loc
                !write(*,*) "dt (yr): ", dtnew(ilevel)*scale_t/(3600*24*365.25), " vol(kpc^3): ", vol_loc*scale_l**ndim
                !if((trel>0.).and.(.not.cosmo)) PoissMean = PoissMean*min((t/trel),1.0)
#if NDIM==2
                PoissMean = PoissMean*4.0*Rad_cloud/3.0
#endif
                ! Compute Poisson realisation
                call poissdev(localseed,PoissMean,nSN)
                if (PoissMean > maxPoissMean) maxPoissMean = PoissMean
!               if (nosn) then
#ifdef SN_INJECT
                if ((firstsn) .and. (x > 0.5*boxlen - 0.75*dx_loc) .and. (y > 0.5*boxlen - 0.75*dx_loc) .and. (x < 0.5*boxlen) .and. (y < 0.5*boxlen) .and. (z > 0.5*boxlen - 0.75*dx_loc) .and. (z < 0.5*boxlen))then
                   nSN = 1
                   firstsn = .false.
                endif
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
                   mSN_loc(nSN_loc)=10.0*2d33/(scale_d*scale_l**3)
                   !if (outputSN) then
                      fileloc=trim(output_dir)//'sn.dat'
                      ilun=140
                      inquire(file=fileloc,exist=file_exist)
                      if(.not.file_exist) then
                         open(ilun, file=fileloc, form='formatted')
                         write(ilun,*)"Time                      x                         y                         z                          Level ProbSN                    rho_sfgas                 ProcID"
                      else
                         open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
                      endif
#if NDIM==2
                      z = 0.0
#endif
                      write(ilun,'(4E26.16,I5,E26.16,E26.16,I5)') t, x, y, z, ilevel, PoissMean, uold(ind_cell(i),imetal+1), myid
                      close(ilun)
                   !endif
                !   do iSN=1,nSN
                !      uold(ind_cell(i),ndim+2) = uold(ind_cell(i),ndim+2) + 10.0*2d33/(scale_d*scale_l**3)*ESN/vol_loc
                !      write(*,*) "SN explosion!"
                !   end do
                end do
              else
                rho_sfr = 0.0
                nSN = 0
              endif
           endif
        enddo
     end do
  end do

#ifdef SNIA_FEEDBACK
  if (nSNIa > 0) then
write(*,*)'min_r2',min_r2(1,1),xSNIa(1,1),xSNIa(1,2),myid
    allocate(min_r2_all(1:2,1:nSNIa))
    call MPI_ALLREDUCE(min_r2,min_r2_all,nSNIa,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,info)
  endif
  do iSN=1,nSNIa
     if (int(min_r2_all(2,iSN)) == myid) then ! The cell containing the SN center is on this processor
       nSN_loc=nSN_loc+1
       xSN_loc(nSN_loc,1)=xSNIa(iSN,1)
       xSN_loc(nSN_loc,2)=xSNIa(iSN,2)
#if NDIM==3
       xSN_loc(nSN_loc,3)=xSNIa(iSN,3)
#endif
       mSN_loc(nSN_loc)=10.0*2d33/(scale_d*scale_l**3)
       !if (outputSN) then
       fileloc=trim(output_dir)//'snIa.dat'
       ilun=140
       inquire(file=fileloc,exist=file_exist)
       if(.not.file_exist) then
          open(ilun, file=fileloc, form='formatted')
          write(ilun,*)"Time                      x                         y                         z ProbSN                    ProcID"
       else
          open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
       endif
#if NDIM==3
       z = xSNIa(iSN,3)
#else
      z = 0.0
#endif
       write(ilun,'(4E26.16,E26.16,I5)') t, xSNIa(iSN,1), xSNIa(iSN,2), z, PoissMeanIa, myid
       close(ilun)
       !endif
      endif
  end do
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
     nsubcycle = nsubcycle(ilevel) ! Assume the same no. of subcycles on every level, will not work otherwise!!!
     weigth = 1d0/(nsubcycle**(ilevel-levelmin))
     sfr_tot(ilevel-levelmin+1) = sfr_tot(ilevel-levelmin+1) + weight*sfr_tot_level
  endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  if (nSN_tot .eq. 0) return

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN explosions=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN(1:nSN_tot,1:ndim))
  allocate(mSN(1:nSN_tot))
  xSN=0.;mSN=0.;
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
    end do
  end do

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:ndim),mSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN,xSN_all,nSN_tot*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN=xSN_all
  mSN=mSN_all
  deallocate(xSN_all,mSN_all)
#endif

  nSN=nSN_tot

  allocate(vol_gas(1:nSN,1:RADCELL_MAX),dq(1:nSN,1:ndim),ekBlast(1:nSN),rSN(1:nSN))
  allocate(indSN(1:nSN))

#ifdef DELAYED_SN
  if(nSN_prev > 0)then
     ! Add SN from previous time step
     ! Compute blast radius
     call subgrid_average_SN(sn_coords,rSN,vol_gas,dq,ekBlast,indSN,nSN_prev)

     ! Modify hydro quantities to account for a Sedov blast wave
     call subgrid_Sedov_blast(sn_coords,mSN,rSN,indSN,vol_gas,dq,ekBlast,nSN_prev)
  endif
#else
  ! Compute blast radius
  call subgrid_average_SN(xSN,rSN,vol_gas,dq,ekBlast,indSN,nSN)

  ! Modify hydro quantities to account for a Sedov blast wave
  call subgrid_Sedov_blast(xSN,mSN,rSN,indSN,vol_gas,dq,ekBlast,nSN)
#endif

#ifdef DELAYED_SN
  if (nSN > 0)then
     if (nSN > NMAX_SN) then
        write(*,*)"WARNING: Too many SN for DELAYED_SN module!!!"
     else
        sn_coords(nSN_prev+1:nSN_prev+nSN) = xSN
        nSN_prev = nSN_prev+nSN
        write(*,*)"nSN_prev",nSN_prev," SNcoords:",sn_coords(nSN_prev+1,1),",",sn_coords(nSN_prev+1,2),",",sn_coords(nSN_prev+1,3)
     endif
  endif
#endif

  deallocate(xSN,mSN,indSN,vol_gas,dq,ekBlast,rSN)

#ifdef SNIA_FEEDBACK
  if (nSNIa > 0)deallocate(xpdf,xSNIa,min_r2,min_r2_all)
#endif


  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

!  write(*,*) "Tot sf: ", tot_sf
!write(*,*) "dens crit (d>", d0, "): ", dpass, " temp crit: ", Tpass, " refine crit: ", gpass


!  !---------------------------------------------------------
!  ! Convert hydro variables back to conservative variables
!  !---------------------------------------------------------
!  ncache=active(ilevel)%ngrid
!  do igrid=1,ncache,nvector
!     ngrid=MIN(nvector,ncache-igrid+1)
!     do i=1,ngrid
!        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
!     end do
!     do ind=1,twotondim
!        iskip=ncoarse+(ind-1)*ngridmax
!        do i=1,ngrid
!           ind_cell(i)=iskip+ind_grid(i)
!        end do
!        do i=1,ngrid
!           d=uold(ind_cell(i),1)
!           u=uold(ind_cell(i),2)
!           v=uold(ind_cell(i),3)
!#if NDIM==3
!           w=uold(ind_cell(i),4)
!#endif
!           !e=uold(ind_cell(i),5)*d
!#ifdef SOLVERmhd
!           bx1=uold(ind_cell(i),6)
!           by1=uold(ind_cell(i),7)
!#if NDIM==3
!           bz1=uold(ind_cell(i),8)
!#endif
!           bx2=uold(ind_cell(i),nvar+1)
!           by2=uold(ind_cell(i),nvar+2)
!#if NDIM==3
!           bz2=uold(ind_cell(i),nvar+3)
!#endif
!           !e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
!#endif
!           !e=e+0.5d0*d*(u**2+v**2+w**2)
!#if NENER>0
!           !do irad=0,nener-1
!           !   e=e+uold(ind_cell(i),inener+irad)
!           !end do
!#endif
!           uold(ind_cell(i),1)=d
!           uold(ind_cell(i),2)=d*u
!           uold(ind_cell(i),3)=d*v
!#if NDIM==3
!           uold(ind_cell(i),4)=d*w
!#endif
!           !uold(ind_cell(i),5)=e
!        end do
!        do ivar=imetal,nvar
!           do i=1,ngrid
!              d=uold(ind_cell(i),1)
!              w=uold(ind_cell(i),ivar)
!              uold(ind_cell(i),ivar)=d*w
!           end do
!        end do
!     end do
!  end do

end subroutine subgrid_sn_feedback
!#endif

!################################################################
!################################################################
!################################################################
!################################################################
subroutine subgrid_average_SN(xSN,rSN,vol_gas,dq,ekBlast,ind_blast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
#if NDIM==2
  use cooling_module, only: twopi
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer,parameter::RADCELL_MAX=10
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,iSN,ind,ix,iy,iz,ngrid,iskip,radcells
  integer::i,nx_loc,igrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,u,v,w,u2,v2,w2,dr_cell,massdiff,mindiff
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc
  integer ,dimension(1:nSN)::ind_blast,SNmaxrad,SNmaxrad_all
  real(dp),dimension(1:nSN)::ekBlast,rSN
  real(dp),dimension(1:nSN,1:RADCELL_MAX)::vol_gas,vol_gas_all,mtot,mtot_all
  real(dp),dimension(1:nSN,1:ndim)::xSN,dq,u2Blast
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN)::ekBlast_all,SNmenc,SNvol
  real(dp),dimension(1:nSN,1:ndim)::dq_all,u2Blast_all
#endif
  logical ,dimension(1:nvector),save::ok

  real(dp)::SN_blast_mass = 60.0

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

  ! Maximum radius of the ejecta
  !rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  !rmax=rmax/scale_l
  rmax=2.0d0*dx_min
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1;rSN=0.0;vol_gas_all=0.0;mtot=0.0;mtot_all=0.0;SNmenc=0.0;SNvol=0.0
  SNmaxrad=RADCELL_MAX

  do iSN=1,nSN
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
                          if(sqrt(dr_SN) .lt. 1.001*dx_min*radcells) then
                             if ((ilevel < nlevelmax) .and. (radcells <= SNmaxrad(iSN))) then
                                if (radcells == 1) then
                                   SNmaxrad(iSN) = 1 ! SN is centered on a coarse cell, this issue will be dealt with in subgrid_sedov_blast
                                else
                                   SNmaxrad(iSN) = radcells-1 ! SN radius must be smaller than this to avoid overlapping coarse cells
                                endif
                             endif
                             mtot(iSN,radcells) = mtot(iSN,radcells) + uold(ind_cell(i),1)*vol_loc
                             vol_gas(iSN,radcells) = vol_gas(iSN,radcells) + vol_loc
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
  call MPI_ALLREDUCE(SNmaxrad,SNmaxrad_all,nSN  ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

  do iSN=1,nSN
     mindiff = 1d10
     do radcells=1,RADCELL_MAX
        massdiff = abs(SN_blast_mass - mtot_all(iSN,radcells)*scale_d*scale_l**3/2d33)
write(*,*)'radcells',radcells,' massdiff',massdiff,' mtot',mtot_all(iSN,radcells)*scale_d*scale_l**3/2d33
        if (massdiff < mindiff) then
           mindiff = massdiff
           rSN(iSN) = radcells*dx_min
           SNmenc(iSN) = mtot_all(iSN,radcells)
           SNvol(iSN) = vol_gas_all(iSN,radcells)
        endif
     enddo
     if (SNmaxrad_all(iSN)*dx_min <= rSN(iSN)) then
        if (myid==1)write(*,*)"WARNING: SN should extend ",rSN(iSN), " kpc but can only extend ",SNmaxrad_all(iSN)*dx_min," kpc without overlapping coarse cells!!!"
        rSN(iSN) = SNmaxrad_all(iSN)*dx_min
        SNmenc(iSN) = mtot_all(iSN,SNmaxrad_all(iSN))
        SNvol(iSN) = vol_gas_all(iSN,SNmaxrad_all(iSN))
     endif

     if(myid==1)write(*,*)"SN rad: ", rSN(iSN), " SN rad (cells): ", rSN(iSN)/dx_min, " SN Mass: ", SNmenc(iSN)*scale_d*scale_l**3/2d33, "blast temperature: ",(1d51*(gamma-1.0)/(SNmenc(iSN)*scale_d*scale_l**3))*(0.6*1.66e-24/1.3806e-16)

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
                    if(dr_SN.lt.rSN(iSN))then
                       uold(ind_cell(i),1) = SNmenc(iSN)/SNvol(iSN)
                    endif
                 endif
              end do
           end do
        end do     ! End loop over grids
     end do    ! End loop over levels
  end do  ! End loop over SNe


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas,vol_gas_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq     ,dq_all     ,nSN*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast,u2Blast_all,nSN*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekBlast,ekBlast_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas=vol_gas_all
  dq     =dq_all
  u2Blast=u2Blast_all
  ekBlast=ekBlast_all
#endif
  do iSN=1,nSN
     if(SNvol(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/SNvol(iSN)
        dq(iSN,2)=dq(iSN,2)/SNvol(iSN)
#if NDIM==3
        dq(iSN,3)=dq(iSN,3)/SNvol(iSN)
#endif
        u2Blast(iSN,1)=u2Blast(iSN,1)/SNvol(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/SNvol(iSN)
#if NDIM==3
        u2Blast(iSN,3)=u2Blast(iSN,3)/SNvol(iSN)
#endif
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
#if NDIM==3
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
#else
        w2=0.0
#endif
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine subgrid_average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subgrid_Sedov_blast(xSN,mSN,rSN,indSN,vol_gas,dq,ekBlast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,vol_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:nSN)::mSN,p_gas,d_gas,d_metal,vol_gas,uSedov,ekBlast,rSN
  real(dp),dimension(1:nSN,1:ndim)::xSN,dq
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

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

  ! Maximum radius of the ejecta
  rmax=2.0d0*dx_min!MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  !rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Minimum star particle mass
!  if(m_star < 0d0)then
!     mstar=n_star/(scale_nH*aexp**3)*vol_min
!  else
!     mstar=m_star*mass_sph
!  endif
!  msne_min=mass_sne_min*2d33/(scale_d*scale_l**3)
!  mstar_max=mass_star_max*2d33/(scale_d*scale_l**3)
  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*2d33))/scale_v**2
  do iSN=1,nSN
!     eta_sn2    = eta_sn
!     if(sf_imf)then
!        if(mSN(iSN).le.mstar_max)then
!           if(mSN(iSN).ge.msne_min) eta_sn2 = eta_ssn
!           if(mSN(iSN).lt.msne_min) eta_sn2 = 0.0
!        endif
!     endif
     if(vol_gas(iSN)>0d0)then
        d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
!        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=mSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*mSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*ESN/ekBlast(iSN))
        endif
     else
        d_gas(iSN)=mSN(iSN)/ekBlast(iSN)
        p_gas(iSN)=mSN(iSN)*ESN/ekBlast(iSN)
!        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/ekBlast(iSN)
     endif
!write(*,*)"SN d: ", d_gas(iSN), " p: ", p_gas(iSN), " uSedov: ", uSedov(iSN), " ekBlast: ", ekBlast(iSN), " vol_gas: ", vol_gas(iSN)
  end do

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
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
#if NDIM==3
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
#else
                    dr_SN=dxx**2+dyy**2
#endif
                    if(dr_SN.lt.rSN(iSN))then!rmax2)then
                       ! Compute the mass density in the cell
                       !uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas(iSN)
                       ! Compute the metal density in the cell
                       !if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_metal(iSN)
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=0.0!uSedov(iSN)*(dxx/rmax-dq(iSN,1))!+vSN(iSN,1)
                       v=0.0!uSedov(iSN)*(dyy/rmax-dq(iSN,2))!+vSN(iSN,2)
!#if NDIM==3
!                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))!+vSN(iSN,3)
!#else
                       w=0.0
!#endif
                       ! Add each momentum component of the blast wave to the gas
                       !uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas(iSN)*u
                       !uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas(iSN)*v
#if NDIM==3
                       !uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas(iSN)*w
#endif
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),ndim+2)=uold(ind_cell(i),ndim+2)+0.5*d_gas(iSN)*(u*u+v*v+w*w)+p_gas(iSN)
!write(*,*)"SN d: ", d_gas(iSN), " vx: ", u, " vy: ", v, " deltaE: ", 0.5*d_gas(iSN)*(u*u+v*v+w*w)+p_gas(iSN)
!write(*,*)"SN rho: ", uold(ind_cell(i),1)," temp: ", (uold(ind_cell(i),ndim+2)*(gamma-1.0)-0.5*(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2)/uold(ind_cell(i),1))*scale_t2/uold(ind_cell(i),1)
                    endif
                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

!  do iSN=1,nSN
!     if(vol_gas(iSN)==0d0)then
!        u=0.0!vSN(iSN,1)
!        v=0.0!vSN(iSN,2)
!#if NDIM==3
!        w=0.0!vSN(iSN,3)
!#endif
!        if(indSN(iSN)>0)then
!           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas(iSN)
!           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas(iSN)*u
!           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas(iSN)*v
!#if NDIM==3
!           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas(iSN)*w
!#endif
!           uold(indSN(iSN),ndim+2)=uold(indSN(iSN),ndim+2)+d_gas(iSN)*0.5*(u*u+v*v+w*w)+p_gas(iSN)
!           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_metal(iSN)
!        endif
!     endif
!  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine subgrid_Sedov_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################

