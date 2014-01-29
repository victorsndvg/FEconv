module module_mphtxt
!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 08/01/2014
!
! PUBLIC PROCEDURES:
!   load_mphtxt: loads a mesh from a MPHTXT format file
!   NOTES: 1) mesh must be composed of a single FE type
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report
use module_convers
use module_manage_mphtxt
use module_mesh
use module_pmh
use module_fe_database_pmh
implicit none


contains

!-----------------------------------------------------------------------
! load_mphtxt: read a MPHTXT file
!-----------------------------------------------------------------------
subroutine load_mphtxt(mphtxtfile, nel, nnod, nver, DI_, LNN, LNV, LNE, LNF, nn, mm, nrc, nra, nrv, z, nsd)
character(len=*),         intent(in)  :: mphtxtfile
integer,                  intent(out) :: nel, nnod, nver, DI_, LNN, LNV, LNE, LNF
integer, allocatable,     intent(out) :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64),allocatable, intent(out) :: z(:,:)
type(mphtxt)                          :: u
type(pmh_mesh)                        :: mphtxt_m
type(mfm_mesh), allocatable           :: m(:)
integer                               :: d, DIM



  !inital settings
  call report_option('level', 'stdout')
  !process MPHTXT file
  print*, 'OPEN MPHTXT FILE:',trim(mphtxtfile)
  call open_mphtxt(u, mphtxtfile)
  print*, 'READING MPHTXT FILE:'
  print*, '-------------------------'
  call read_mphtxt(u, m, mphtxt_m, DIM)
!  call counter_clockwise_orientation(mphtxt_m)
  call mphtxt2mfm(mphtxt_m, m)


!  print*, 'DIM', DIM
  !allocate mesh(es)
!  allocate(m(DIM))

!  forall (d = 1:DIM) m(d)%DIM = d
  print*, 6
  !return msh variables
  nel  = m(DIM)%nl
  nnod = m(DIM)%nd
  nver = m(DIM)%nv
  DI_  = m(DIM)%DIM
  LNN  = m(DIM)%LNN
  LNV  = m(DIM)%LNV
  LNE  = m(DIM)%LNE
  LNF  = m(DIM)%LNF
  print*, 7
  if (allocated(m(DIM)%id)) then; allocate( nn(LNN,nel));  nn(1:LNN,1:nel)  = m(DIM)%id; end if
  print*, 8
  if (allocated(m(DIM)%iv)) then; allocate( mm(LNV,nel));  mm(1:LNV,1:nel)  = m(DIM)%iv; end if
  print*, 9
  if (allocated(m(DIM)%rf)) then; allocate(nrc(LNF,nel)); nrc(1:LNF,1:nel)  = m(DIM)%rf; end if
  print*, 10
  if (allocated(m(DIM)%re)) then; allocate(nra(LNE,nel)); nra(1:LNE,1:nel)  = m(DIM)%re; end if
  print*, 11
  if (allocated(m(DIM)%rv)) then; allocate(nrv(LNV,nel)); nrv(1:LNV,1:nel)  = m(DIM)%rv; end if
  print*, 12
  if (allocated(m(DIM)%xv)) then; allocate(  z(DIM,nver));  z=0.0_real64; z(1:DIM,1:nver) = m(DIM)%xv; end if
  print*, 13
  if (allocated(m(DIM)%rl)) then; allocate(    nsd(nel));       nsd(1:nel)  = m(DIM)%rl; end if
  print*, 14


  !deallocate mesh
  do d = 1, DIM
    call dealloc_mesh(m(d))
  enddo
  deallocate(m)

  print*, 15

end subroutine

!-----------------------------------------------------------------------
! mphtxt2mfm: build a MFM mesh from a MPHTXT mesh structure
!-----------------------------------------------------------------------
subroutine mphtxt2mfm(mphtxt_m, m)
  type(pmh_mesh),                          intent(in)    :: mphtxt_m   ! MPHTXT mesh
  type(mfm_mesh), dimension(:),allocatable,intent(inout) :: m          ! MFM mesh
  integer                                                :: maxfetype  ! Higher order FE type
  integer                                                :: maxfelnv    ! Max NN FE 
  integer                                                :: sdim, numpoints, numvertex, numedges, numsurfels, numvolels, numel
  integer                                                :: i, j, aux, iniobj, endobj, inietype, endetype, etype

  if (allocated(m)) deallocate(m)
  sdim = 0
  maxfetype  = 0
  maxfelnv   = 0
  numpoints  = 0
  numvertex  = 0
  numedges   = 0
  numsurfels = 0
  numvolels  = 0
  numel      = 0

  ! loop in objects
  do i = 1, size(mphtxt_m%pc,1)
    if (sdim < mphtxt_m%pc(i)%dim) sdim = mphtxt_m%pc(i)%dim
    print*, 'dimension: sdim:', sdim, 'mphtxt%dim:',mphtxt_m%pc(i)%dim
    numpoints = numpoints + mphtxt_m%pc(i)%nnod
    do j = 1, size(mphtxt_m%pc(i)%el,1)
      print*, '::::::::::::::::::::::::::::::::::'
      print*, 'ELEMT TYPE', mphtxt_m%pc(i)%el(j)%type
      print*, '::::::::::::::::::::::::::::::::::'
      if (maxfelnv < FEDB(mphtxt_m%pc(i)%el(j)%type)%lnv) then 
        maxfetype = mphtxt_m%pc(i)%el(j)%type
        maxfelnv = FEDB(mphtxt_m%pc(i)%el(j)%type)%lnv
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Only tetra, triangle, edge and vertex (P1 and P2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      etype = mphtxt_m%pc(i)%el(j)%type
      if (etype == NODE2 .or. etype == NODE3) then       ! Vertex
        numvertex = numvertex + mphtxt_m%pc(i)%el(j)%nel
      elseif (etype == ED2_P1 .or. etype == ED3_P1) then ! Edge P1
        numedges = numedges + mphtxt_m%pc(i)%el(j)%nel
      elseif (etype == TR2_P1 .or. etype == TR3_P1 .or. etype == QU2_P1 .or. etype == QU3_P1) then ! Triangle P1
        numsurfels = numsurfels + mphtxt_m%pc(i)%el(j)%nel;
      elseif (etype == TET_P1 .or. etype == HEX_P1) then                      ! Tetrahedra P1
        numvolels = numvolels + mphtxt_m%pc(i)%el(j)%nel
      else
        call error("Finite element type not allowed in Modulef")
      endif
    enddo
  enddo

  if (numvolels /= 0) then
    numel = numvolels
  elseif (numsurfels/= 0) then
    numel = numsurfels
  elseif (numedges/= 0) then
    numel = numedges
  elseif (numpoints/= 0) then
    numel = numpoints
  else
    call error('mphtxt/mphtxt2mfm: Wrong number of elements')
  endif

  ! MFM mesh initialization
  allocate(m(sdim))
  m(sdim)%DIM = sdim                               ! Space dimension
  m(sdim)%LNN = FEDB(maxfetype)%lnn              ! Local number of nodes
  m(sdim)%LNV = FEDB(maxfetype)%lnv              ! Local number of vertices
  m(sdim)%LNE = FEDB(maxfetype)%lne              ! Local number of elements
  m(sdim)%LNF = FEDB(maxfetype)%lnf              ! Local number of faces

  allocate(m(sdim)%xv(sdim,numpoints))     ! z
  allocate(m(sdim)%id(m(sdim)%LNN, numel)) ! nn
  allocate(m(sdim)%iv(m(sdim)%LNV, numel)) ! mm
  allocate(m(sdim)%rv(m(sdim)%LNV, numel)) ! nrv
  allocate(m(sdim)%re(m(sdim)%LNE, numel)) ! nra
  allocate(m(sdim)%rf(m(sdim)%LNF, numel)) ! nrc
  allocate(m(sdim)%rl(numel))              ! nsd

  print*, 'numvolels:', numvolels
  print*, 'numsurfels:',numsurfels
  print*, 'numedges:', numedges
  print*, 'numpoints:', numpoints

  m(sdim)%nl = numel !global number of elements
  m(sdim)%nd = numpoints !global number of nodes
  m(sdim)%nv = numpoints !global number of vertices
  m(sdim)%rv(:,:) = 0
  m(sdim)%re(:,:) = 0
  m(sdim)%rf(:,:) = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*, 'size(m):   ',size(m)
  print*, 'size(z):   ',size(m(sdim)%xv,1), size(m(sdim)%xv,2)
  print*, 'size(nn):  ',size(m(sdim)%id,1), size(m(sdim)%id,2)
  print*, 'size(mm):  ',size(m(sdim)%iv,1), size(m(sdim)%iv,2)
  print*, 'size(nrv): ',size(m(sdim)%rv,1), size(m(sdim)%rv,2)
  print*, 'size(nra): ',size(m(sdim)%re,1), size(m(sdim)%re,2)
  print*, 'size(nrc): ',size(m(sdim)%rf,1), size(m(sdim)%rf,2)
  print*, 'size(nsd): ',size(m(sdim)%rl,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*, 1
  call mphtxt_build_nn_nsd(mphtxt_m, maxfetype, sdim, m)
  print*, 2
  if (numvolels /= 0) call mphtxt_build_nrc(mphtxt_m, maxfetype, sdim, m)
  print*, 3
  if (numsurfels /= 0) call mphtxt_build_nra(mphtxt_m, maxfetype, sdim, m)
  print*, 4
  call mphtxt_build_nrv(mphtxt_m, maxfetype, sdim, m)
  print*, 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  do i = 1, numel
!!    print*, m(sdim)%iv(:,i),'----', m(sdim)%rl(i)
!!  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine


!!-----------------------------------------------------------------------
!! mphtxt_build_nn: build nodes connectivity and groups (nsd)
!!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   TEMPORAL: Construye tambien MM, SOLO en Lagrange P1 Coincide!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mphtxt_build_nn_nsd(mphtxt_m, maxfetype, sdim, m)
  type(pmh_mesh),                             intent(in)    :: mphtxt_m   ! MPHTXT mesh
  type(mfm_mesh), dimension(:), allocatable, intent(inout) :: m          ! MFM mesh
  integer,                                intent(in)    :: maxfetype  ! Higher order FE type 
  integer,                                intent(in)    :: sdim       ! Space dimension
  integer                                               :: i, j, aux, iniobj, endobj, inietype, endetype

  iniobj = 1
  endobj = 1
  inietype = 1
  endetype = 1

  do i = 1, size(mphtxt_m%pc,1)
    endobj = iniobj + mphtxt_m%pc(i)%nnod - 1
    m(sdim)%xv(:,iniobj:endobj) = mphtxt_m%pc(i)%znod(:,:)
    iniobj = endobj + 1
    do j = 1, size(mphtxt_m%pc(i)%el,1)
      if (mphtxt_m%pc(i)%el(j)%type == maxfetype) then
        endetype = inietype + mphtxt_m%pc(i)%el(j)%nel -1
        m(sdim)%id(:,inietype:endetype) = mphtxt_m%pc(i)%el(j)%nn ! nn
        m(sdim)%iv(:,inietype:endetype) = mphtxt_m%pc(i)%el(j)%nn ! mm: wrong!!!!!!!
        m(sdim)%rl(inietype:endetype) = mphtxt_m%pc(i)%el(j)%ref(:)
        inietype = endetype + 1 
      endif
    enddo
  enddo

end subroutine

!!-----------------------------------------------------------------------
!! mphtxt_build_nrc: build nrc. Must be called after build mm
!!-----------------------------------------------------------------------
subroutine mphtxt_build_nrc(mphtxt_m, maxfetype, sdim, m)
  type(pmh_mesh),                             intent(in)    :: mphtxt_m   ! MPHTXT mesh
  type(mfm_mesh), dimension(:), allocatable, intent(inout) :: m          ! MFM mesh
  integer,                                intent(in)    :: maxfetype  ! Higher order FE type 
  integer,                                intent(in)    :: sdim       ! Space dimension
  integer                                               :: i, j, ii, jj, kk, indx
  integer, dimension(3)                                 :: pos
  logical                                               :: founded

  do i = 1, size(mphtxt_m%pc,1)                              ! loop over all mphtxt objects
    do j = 1, size(mphtxt_m%pc(i)%el,1)                      ! loop over all element types
      if (mphtxt_m%pc(i)%el(j)%type == TR2_P1 .or. mphtxt_m%pc(i)%el(j)%type == TR3_P1) then   ! Check element vertices (3 = triangles)
	print*, sdim, size(m(sdim)%iv,1), size(m(sdim)%iv,2), size(mphtxt_m%pc(i)%el(j)%nn,1),size(mphtxt_m%pc(i)%el(j)%nn,2)
        do ii = 1, size(m(sdim)%iv,2)                        ! loop over mesh connectivity array (mm)
          do jj = 1, mphtxt_m%pc(i)%el(j)%nel                ! loop over all elements in a type
            founded = .true.
            do kk = 1, size(mphtxt_m%pc(i)%el(j)%nn,1)       ! loop over all nodes in a element
              if (.not. any(m(sdim)%iv(:,ii) == (mphtxt_m%pc(i)%el(j)%nn(kk, jj)),1)) then
                founded = .false.
                exit
              endif
            enddo
            if (founded) then
              call search_pos(m(sdim)%iv(:,ii), mphtxt_m%pc(i)%el(j)%nn(:, jj), pos, size(m(sdim)%iv,1), size(mphtxt_m%pc(i)%el(j)%nn,1))
              call sort_bubble(pos, size(mphtxt_m%pc(i)%el(j)%nn,1))
              indx = mphtxt_m%pc(i)%el(j)%ref(jj) 
              if (pos(1) == 1 .and. pos(2) == 2 .and. pos(3) == 3) then; m(sdim)%rf(1, ii) = indx;
              elseif (pos(1) == 1 .and. pos(2) == 3 .and. pos(3) == 4) then; m(sdim)%rf(2, ii) = indx;
              elseif (pos(1) == 1 .and. pos(2) == 2 .and. pos(3) == 4) then; m(sdim)%rf(3, ii) = indx;
              elseif (pos(1) == 2 .and. pos(2) == 3 .and. pos(3) == 4) then; m(sdim)%rf(4, ii) = indx;
              endif
            endif
          enddo
        enddo        
      endif
    enddo
  enddo

end subroutine

!!-----------------------------------------------------------------------
!! mphtxt_build_nra: build nra. Must be called after build mm
!!-----------------------------------------------------------------------
subroutine mphtxt_build_nra(mphtxt_m, maxfetype, sdim, m)
  type(pmh_mesh),                            intent(in)    :: mphtxt_m   ! MPHTXT mesh
  type(mfm_mesh), dimension(:), allocatable, intent(inout) :: m          ! MFM mesh
  integer,                                   intent(in)    :: maxfetype  ! Higher order FE type 
  integer,                                   intent(in)    :: sdim       ! Space dimension
  integer                                                  :: i, j, ii, jj, kk, indx
  integer, dimension(2)                                    :: pos
  logical                                                  :: founded

  do i = 1, size(mphtxt_m%pc,1)                               ! loop over all mphtxt objects
    do j = 1, size(mphtxt_m%pc(i)%el,1)                       ! loop over all element types
      if (mphtxt_m%pc(i)%el(j)%type == ED2_P1 .or. mphtxt_m%pc(i)%el(j)%type == ED3_P1 ) then      ! Check element dimension (triangles)
        do ii = 1, size(m(sdim)%iv,2)                         ! loop over mesh connectivity array (mm)
          do jj = 1, mphtxt_m%pc(i)%el(j)%nel                 ! loop over all elements in a type
            founded = .true.
            do kk = 1, size(mphtxt_m%pc(i)%el(j)%nn,1)        ! loop over all nodes in a element
              if (.not. any(m(sdim)%iv(:,ii) == (mphtxt_m%pc(i)%el(j)%nn(kk, jj)),1)) then
                founded = .false.
                exit
              endif
            enddo
            if (founded) then
              call search_pos(m(sdim)%iv(:,ii), mphtxt_m%pc(i)%el(j)%nn(:, jj), pos, size(m(sdim)%iv,1), size(mphtxt_m%pc(i)%el(j)%nn,1))
              call sort_bubble(pos, size(mphtxt_m%pc(i)%el(j)%nn,1))
              indx = mphtxt_m%pc(i)%el(j)%ref(jj)
              if (maxfetype == TR2_P1 .or. maxfetype == TR3_P1 .or. maxfetype == TET_P1) then
                if (pos(1) == 1 .and. pos(2) == 2) then; m(sdim)%re(1, ii) = indx;
                elseif (pos(1) == 2 .and. pos(2) == 3) then; m(sdim)%re(2, ii) = indx;
                elseif (pos(1) == 1 .and. pos(2) == 3) then; m(sdim)%re(3, ii) = indx;
                elseif (pos(1) == 1 .and. pos(2) == 4) then; m(sdim)%re(4, ii) = indx;
                elseif (pos(1) == 2 .and. pos(2) == 4) then; m(sdim)%re(5, ii) = indx;
                elseif (pos(1) == 3 .and. pos(2) == 4) then; m(sdim)%re(6, ii) = indx;
                endif
              elseif (maxfetype == QU2_P1 .or. maxfetype == QU3_P1) then
                if (pos(1) == 1 .and. pos(2) == 2) then; m(sdim)%re(1, ii) = indx;
                elseif (pos(1) == 2 .and. pos(2) == 3) then; m(sdim)%re(2, ii) = indx;
                elseif (pos(1) == 3 .and. pos(2) == 4) then; m(sdim)%re(3, ii) = indx;
                elseif (pos(1) == 1 .and. pos(2) == 4) then; m(sdim)%re(4, ii) = indx;
                endif
              else 
                call error('module_mphtxt/mphtxt_build_nra# Element type not supported yet')	
              endif
            endif
          enddo
        enddo        
      endif
    enddo
  enddo

end subroutine

!!-----------------------------------------------------------------------
!! mphtxt_build_nrv: build nrv. Must be called after build mm
!!-----------------------------------------------------------------------
subroutine mphtxt_build_nrv(mphtxt_m, maxfetype, sdim, m)
  type(pmh_mesh),                             intent(in)    :: mphtxt_m   ! MPHTXT mesh
  type(mfm_mesh), dimension(:), allocatable, intent(inout) :: m          ! MFM mesh
  integer,                                intent(in)    :: maxfetype  ! Higher order FE type 
  integer,                                intent(in)    :: sdim       ! Space dimension
  integer                                               :: i, j, ii, jj, kk, indx
  integer, dimension(1)                                 :: pos
  logical                                               :: founded

  do i = 1, size(mphtxt_m%pc,1)                                            ! loop over all mphtxt objects
    do j = 1, size(mphtxt_m%pc(i)%el,1)                              ! loop over all element types
      if (mphtxt_m%pc(i)%el(j)%type == NODE2 .or. mphtxt_m%pc(i)%el(j)%type == NODE3) then ! Check element dimension (triangles)
        do ii = 1, size(m(sdim)%iv,2)                                  ! loop over mesh connectivity array (mm)
          do jj = 1, mphtxt_m%pc(i)%el(j)%nel           ! loop over all elements in a type
            founded = .true.
            do kk = 1, size(mphtxt_m%pc(i)%el(j)%nn,1)      ! loop over all nodes in a element
              if (.not. any(m(sdim)%iv(:,ii) == (mphtxt_m%pc(i)%el(j)%nn(kk, jj)),1)) then
                founded = .false.
                exit
              endif
            enddo
            if (founded) then
              call search_pos(m(sdim)%iv(:,ii), mphtxt_m%pc(i)%el(j)%nn(:, jj), pos, size(m(sdim)%iv,1), size(mphtxt_m%pc(i)%el(j)%nn,1))
              indx = mphtxt_m%pc(i)%el(j)%ref(jj)
              if (pos(1) == 1) then; m(sdim)%rv(1, ii) = indx;
              elseif (pos(1) == 2) then; m(sdim)%rv(2, ii) = indx;
              elseif (pos(1) == 3) then; m(sdim)%rv(3, ii) = indx;
              elseif (pos(1) == 4) then; m(sdim)%rv(4, ii) = indx;
              endif
            endif
          enddo
        enddo        
      endif
    enddo
  enddo

end subroutine

!-----------------------------------------------------------------------
! search_pos: Search values array positions in origin array
!-----------------------------------------------------------------------
subroutine search_pos(origin, values, pos, adim, vdim)
  integer,                                intent(in)    :: adim, vdim     ! Array and values dimensions
  integer, dimension(adim),               intent(in)    :: origin         ! Array where to find
  integer, dimension(vdim),               intent(in)    :: values         ! Array of values to find
  integer, dimension(vdim),               intent(out)   :: pos            ! Array with the positions of values
  integer                                               :: i, j

  if (size(values,1) /= size(pos,1)) call error('search_pos/ Array sizes must agree #' //string(size(values,1))//'-'//string(size(pos,1)))
  pos = 0
  do i = 1, vdim
    do j = 1, adim
      if (values(i) == origin(j)) pos(i) = j
    enddo
  enddo

end subroutine

!-----------------------------------------------------------------------
! sort_buble: sorts an array of dimension adim
!-----------------------------------------------------------------------
subroutine sort_bubble (array, adim)
  integer,                  intent(in)    :: adim
  integer, dimension(adim), intent(inout) :: array
  integer                                 :: i = 0
  integer                                 :: j = 0
  integer                                 :: temp

  do i = adim -1, 0, -1
    do j = 1, i, 1
      if (array(j) > array(j+1))  then
        temp = array(j+1)
        array(j+1) = array(j)
        array(j) = temp
      endif
    enddo
  enddo

end subroutine

end module
