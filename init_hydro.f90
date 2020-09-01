subroutine init_hydro
  use amr_commons
  use hydro_commons
  use amr_parameters
#ifdef RT
  use rt_parameters,only: convert_birth_times
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=256)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1108
#if NENER>0
  integer::irad
#endif
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  real(dp)::rad,vel,x1c,x2c,x3c,dx,scale
  integer::ix,iy,iz,nx_loc

  if(verbose)write(*,*)'Entering init_hydro'

  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0
  if(momentum_feedback)then
     allocate(pstarold(1:ncell))
     allocate(pstarnew(1:ncell))
     pstarold=0.0d0; pstarnew=0.0d0
  endif
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+103
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc=trim(output_dir)//'output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc=trim(output_dir)//'output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     endif



     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(myid==1)then
        write(*,*)'Restart - Non-thermal pressure / Passive scalar mapping'
        write(*,'(A50)')"__________________________________________________"
        do i=1,nvar2-(ndim+2)
            if(remap_pscalar(i).gt.0) then
               write(*,'(A,I3,A,I3)') ' Restart var',i+ndim+2,' loaded in var',remap_pscalar(i)
            else if(remap_pscalar(i).gt.-1)then
               write(*,'(A,I3,A)') ' Restart var',i+ndim+2,' read but not loaded'
            else
               write(*,'(A,I3,A)') ' Restart var',i+ndim+2,' not read'
            endif
        enddo
        write(*,'(A50)')"__________________________________________________"
     endif
#ifdef RT
     if((neq_chem.or.rt).and.nvar2.lt.nvar)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        if(myid==1) write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to zero'
     end if
     if((neq_chem.or.rt).and.nvar2.gt.nvar)then ! Not OK to drop variables
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found   =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        call clean_stop
     end if
#endif
     x1c = x1_c*boxlen
     x2c = x2_c*boxlen
     x3c = x3_c*boxlen
!     open(120,file=trim(output_dir)//'coords.dat'//TRIM(nchar),status="new")

     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)

     do ilevel=1,nlevelmax2
        ! Mesh size at level ilevel
        dx=0.5D0**ilevel
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 ! Set position of cell centers relative to grid center
                 iz=(ind-1)/4
                 iy=(ind-1-4*iz)/2
                 ix=(ind-1-2*iy-4*iz)
                 if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
                 if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
                 if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx

                 ! Read density and velocities --> density and momenta
                 do ivar=1,ndim+1
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       do i=1,ncache
                          vel = 0
                          rad = sqrt(((xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale - x1c)**2 + ((xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale - x2c)**2 + ((xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale - x3c)**2)
!                          if (ivar==3)write(120,*)rad,xg(ind_grid(i),1),xc(ind,1),x1c,xg(ind_grid(i),2),xc(ind,2),x2c,xg(ind_grid(i),3),xc(ind,3),x3c
                          if ((ivar==3) .and. (rad > rad_wind) .and. (uold(ind_grid(i)+iskip,1) < rhomax_wind)) vel=vel_wind
                          uold(ind_grid(i)+iskip,ivar)=(xx(i)+vel)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 end do

#if NENER>0
                 ! Read non-thermal pressures --> non-thermal energies
                 do ivar=ndim+3,ndim+2+nener
                    if(remap_pscalar(ivar-ndim-2).gt.-1) read(ilun)xx
                    do i=1,ncache
                       if(remap_pscalar(ivar-ndim-2).gt.0) then
                          uold(ind_grid(i)+iskip,remap_pscalar(ivar-ndim-2))=xx(i)/(gamma_rad(ivar-ndim-2)-1d0)
                       else if(remap_pscalar(ivar-ndim-2).lt.0) then
                          uold(ind_grid(i)+iskip,abs(remap_pscalar(ivar-ndim-2)))=0d0
                       endif
                    end do
                 end do
#endif
                 ! Read thermal pressure --> total fluid energy
                 read(ilun)xx
                 do i=1,ncache
                    xx(i)=xx(i)/(gamma-1d0)
                    if (uold(ind_grid(i)+iskip,1)>0.)then
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NENER>0
                    do irad=1,nener
                       xx(i)=xx(i)+uold(ind_grid(i)+iskip,ndim+2+irad)
                    end do
#endif
                 else
                    xx(i)=0.
                 end if
                    uold(ind_grid(i)+iskip,ndim+2)=xx(i)
                 end do
#if NVAR>NDIM+2+NENER
                 ! Read passive scalars
                 do ivar=ndim+3+nener,max(nvar2,nvar)
                    if(remap_pscalar(ivar-ndim-2).gt.-1) read(ilun)xx
                    if(ivar.gt.nvar)then
                       continue
                    endif
                    do i=1,ncache
                       if(remap_pscalar(ivar-ndim-2).gt.0)then
                          uold(ind_grid(i)+iskip,remap_pscalar(ivar-ndim-2))=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       else if(remap_pscalar(ivar-ndim-2).lt.0) then
                          uold(ind_grid(i)+iskip,abs(remap_pscalar(ivar-ndim-2)))=0d0
                       endif
                    end do
                 end do
#endif
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)
!     close(120)
!     call clean_stop

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif



#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro




