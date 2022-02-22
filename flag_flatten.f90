subroutine flag_flatten(ilevel)
  ! Flag cells with velocity+sound speed (in each direction) above threshold and their immediate neighbors (in that direction)
  ! The HLL solver and, optionally, the minmod slope limiter will then be used in such cells rather than e.g. HLLC and moncen to add extra artificial diffusion
  ! This is to help smooth out such "problem cells" to alleviate issues with extremely short time steps and negative densities and pressures
  ! This method is akin to the PLUTO code's SHOCK_FLATTENING implementation, although less sophisticated
  ! Flagging is done by changing metallicity to save memory artificialy boosting it by a factor of 1e6
  ! Assumes 'real' metallicities are always within the range 1e-5 < Z/Z_solar < 10
  ! METALLICITIES MUST THEN BE FIXED DURING THE 'UNSPLIT' SUBROUTINE IN UMUSCL.F90!

  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,ncache,i,j,iskip
  integer::igrid,ind,idim,ngrid,ivar
  integer::nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  integer,dimension(1:nvector,1:twondim),save::indn
  real(dp)::vel,cs,prs

  ! Loop over active grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring offsets
     call getnborgrids(ind_grid,igridn,ngrid)

     ! Loop over cells
     do ind=1,twotondim

        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather neighboring cells
        call getnborcells(igridn,ind,indn,ngrid)

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do j=1,twondim
           do i=1,ngrid
              if(indn(i,j)==0)then
                 indn(i,j)=nbor(ind_grid(i),j)
              end if
           end do
        end do

        ! Loop over dimensions
!        do idim=1,ndim
           ! Gather hydro variables
        do i=1,ngrid
              vel = sqrt((uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2)/max(uold(ind_cell(i),1),smallr)**2)
              prs=uold(ind_cell(i),ndim+2)-0.5d0*max(uold(ind_cell(i),1),smallr)*vel**2
#if NENER>0
              do irad=0,nener-1
                 prs=prs-uold(ind_cell(i),inener+irad)
              end do
#endif
              prs = prs*(gamma-1.0)
              cs = sqrt(gamma*prs/max(uold(ind_cell(i),1),smallr))
              if (cs+vel > sound_speed_thresh) then ! Flag for flattening
                 write(*,*)'Flagged cell for flattening with c_s=',cs,'vel=',vel,'rho=',uold(ind_cell(i),1),'prs=',prs,'T=',0.6*prs/(uold(ind_cell(i),1)*7.75d-5)
                 if (uold(ind_cell(i),imetal)/uold(ind_cell(i),1) < 10.0)uold(ind_cell(i),imetal) = uold(ind_cell(i),imetal)*1d6 ! If not flagged already
                 do idim=1,ndim
                    if (uold(indn(i,2*idim-1),imetal)/uold(indn(i,2*idim-1),1) < 10.0)uold(indn(i,2*idim-1),imetal) = uold(indn(i,2*idim-1),imetal)*1d6
                    if (uold(indn(i,2*idim  ),imetal)/uold(indn(i,2*idim  ),1) < 10.0)uold(indn(i,2*idim  ),imetal) = uold(indn(i,2*idim  ),imetal)*1d6
                 enddo
              endif
           end do
!      	end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine flag_flatten
