module module_read_mphtxt
!-----------------------------------------------------------------------
! Module for mphtxt file read
! Last update: 08/01/2014
!-----------------------------------------------------------------------
use module_COMPILER_DEPENDANT, only: real64
use module_os_dependant, only: maxpath
use module_report, only:error
use module_convers
use module_ALLOC
use module_mesh
use module_pmh
contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! read: read mphtxt file header
!-----------------------------------------------------------------------

subroutine read_mphtxt_header(iu, mphtxt_m)

  integer,                          intent(in)    :: iu             ! Unit number for mphtxtfile
  type(pmh_mesh),                   intent(inout) :: mphtxt_m
  integer,dimension(2)                            :: version        ! Version of mphtxtfile
  character(len=MAXPATH)                          :: line
  integer                                         :: aux, i, j
  integer                                         :: ntags, ntypes
  character(len=MAXPATH),dimension(:),allocatable :: tags, types    ! Tags and typess

    ! Variables initialization
    if (allocated(tags)) deallocate(tags)
    if (allocated(mphtxt_m%pc)) deallocate(mphtxt_m%pc)
    version = (/-1,-1/)
    ntags = -1
    ntypes = -1
    i = 1
    j = 1

    do 
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios /= 0) call error('mphtxt_file/header, #'//trim(string(ios)))

      line = trim(line,'#')

      if (len_trim(line) /= 0) then ! Discards empty lines
        ! File version
        if (version(1) == -1 .and. version(2) == -1 .and. word_count(line,' ') == 2) then
          read(line,*) version(1),version(2)
          if (version(1) /= 0 .or. version(2) /= 1) call error('mphtxt_file/header, 0.1 is the only supported version #'//&
          &string(version))  
        ! Number of tags
        elseif (ntags == -1 .and. word_count(line,' ') == 1) then
          read(line,*) ntags
          allocate(tags(ntags))
        ! Tags
        elseif (allocated(tags) .and. (i <= ntags) .and. (word_count(line,' ') == 2)) then
          read(line,*) aux,tags(i)
          i = i+1
          if (i > ntags) cycle ! All tags already readed
        ! Number of types
        elseif (ntypes == -1 .and. word_count(line,' ') == 1) then
          read(line,*) ntypes
          allocate(types(ntypes))
          allocate(mphtxt_m%pc(ntypes))
        ! Types and objects
        elseif (allocated(types) .and. (j <= ntypes) .and. (word_count(line,' ') == 2)) then
          read(line,*) aux,types(j)
          j = j+1
          if (j > ntypes) exit ! All types already readed. Header readed.
        else
          exit
        endif
      endif
    enddo


end subroutine


!-----------------------------------------------------------------------
! read: read a object from mphtxt file
!-----------------------------------------------------------------------
subroutine read_mphtxt_object(iu, mphtxt_o)

  integer,                            intent(in)    :: iu            ! Unit number for mphtxtfile
  type(piece),                        intent(inout) :: mphtxt_o      ! mphtxt mesh
  integer, dimension(3)                             :: serializable  ! Serializable option
  character(len=MAXPATH)                            :: obj_class     ! Object class. Must by 'Mesh'.
  integer                                           :: version       ! Object version
  character(len=MAXPATH)                            :: line
  integer                                           :: aux, i
  integer                                           :: netypes       ! Number of element types
  integer                                           :: offset        ! lowest number of node

    ! Variables initialization
    if (allocated(mphtxt_o%z)) deallocate(mphtxt_o%z)
!    if (allocated(mphtxt_o%etypes)) deallocate(mphtxt_o%etypes)
    mphtxt_o%dim      = -1
    mphtxt_o%nnod     = -1
    offset            = -1
    netypes           = -1 
    serializable      = (/-1,-1,-1/)
    obj_class         = ''
    version           = -1
    i = 1

    ! Read object header
    do 
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios /= 0) call error('mphtxt_file/object, #'//trim(string(ios)))

      line = trim(line,'#')

      if (len_trim(line) /= 0) then  ! Discards empty lines
        ! Serializable object
        if (serializable(1) == -1 .and. serializable(2) == -1 .and. serializable(3) == -1 .and. word_count(line,' ') == 3) then
          read(line,*) serializable(1), serializable(2), serializable(3)
          if (.not. (serializable(1) == 0 .or. serializable(1) == 1)) call error('mphtxt_file/object, Only versions 0 or 1 &
          &supported #'//string(serializable))
          if (serializable(3) /= 1) call error('mphtxt_file/object, Not serializable object #'//string(serializable))
        ! Object class
        elseif (obj_class == '' .and. word_count(line,' ') == 2) then
          read(line,*) aux, obj_class
          if (obj_class /= 'Mesh') call error('mphtxt_file/object, Only mesh object allowed #'//trim(obj_class))
        ! Object version
        elseif (version == -1 .and. (word_count(line,' ') == 1)) then
          read(line,*) version
        ! Space dimension
        elseif ((mphtxt_o%dim == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) mphtxt_o%dim
        ! Number of points
        elseif ((mphtxt_o%nnod == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) mphtxt_o%nnod
        ! Lowest mesh point index
        elseif ((offset == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) offset
          allocate(mphtxt_o%z(mphtxt_o%dim,mphtxt_o%nnod))
        ! Point coordinates
        elseif (allocated(mphtxt_o%z) .and. (i <= mphtxt_o%nnod) .and. (word_count(line,' ') == mphtxt_o%dim)) then
             read(line,*) (mphtxt_o%z(k,i), k=1,mphtxt_o%dim)
             i = i+1
             if (i > mphtxt_o%nnod) cycle ! All coords already readed
        elseif ((netypes == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) netypes
          allocate(mphtxt_o%el(netypes))
          exit ! Number of elements already readed. Object header readed.
        else
          exit
        endif
      endif
    enddo

    ! Read object element types
    do i = 1, netypes
      call read_mphtxt_etype(iu, mphtxt_o%el(i), offset)
    enddo

end subroutine


!type elgroup
!  integer              :: type = 0 !element type (one of those defined in module_eltype)
!  integer              :: nel  = 0 !total number of elements
!  integer, allocatable :: nn(:)    !global numbering of nodes
!  integer, allocatable :: mm(:)    !global numbering of vertices
!  integer, allocatable :: ref(:)   !reference numbering
!end type

!-----------------------------------------------------------------------
! read: read a element type from mphtxt file
!-----------------------------------------------------------------------
subroutine read_mphtxt_etype(iu, mphtxt_t, offset)

  integer,                            intent(in)    :: iu            ! Unit number for mphtxtfile
  type(elgroup),                      intent(inout) :: mphtxt_t      ! mphtxt mesh
  integer,                            intent(in)    :: offset        ! Lowest number of node
  character(len=MAXPATH)                            :: fetype_name   ! Object class. Must by 'Mesh'.
  character(len=MAXPATH)                            :: line
  integer                                           :: aux, i, j,k ,l
  integer                                           :: local_nparam, nparam, local_nnodes, nupdownpairs, nindices

    ! Variables initialization
    if (allocated(mphtxt_t%nn)) deallocate(mphtxt_t%nn)
    if (allocated(mphtxt_t%ref)) deallocate(mphtxt_t%ref)
    fetype_name           = ''
    mphtxt_t%nel          = -1
    local_nnodes          = -1
    nindices              = -1
    local_nparam          = -1
    nparam                = -1
    nupdownpairs          = -1
    i = 1
    j = 1
    k = 1
    l = 1

    ! Read object header
    do 
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios /= 0) call error('mphtxt_file/object/etype, #'//trim(string(ios)))

      line = trim(line,'#')

      if (len_trim(line) /= 0) then  ! Discards empty lines
        ! FE type
        if (len_trim(fetype_name) == 0 .and. word_count(line,' ') == 2) then
          read(line,*) aux, fetype_name
          mphtxt_t%type = mphtxt_get_type(fetype_name)
        ! Local number of nodes per element
        elseif (local_nnodes == -1 .and. word_count(line,' ') == 1) then
          read(line,*) local_nnodes
        ! Number of elements
        elseif (mphtxt_t%nel == -1 .and. word_count(line,' ') == 1) then
          read(line,*) mphtxt_t%nel
          allocate(mphtxt_t%nn(local_nnodes, mphtxt_t%nel))
        ! Elements
        elseif (allocated(mphtxt_t%nn) .and. (i <= mphtxt_t%nel).and. word_count(line,' ') == local_nnodes) then
          read(line,*) (mphtxt_t%nn(m,i), m=1,local_nnodes)
          mphtxt_t%nn(:,i) = mphtxt_t%nn(:,i) - offset + 1
          call mphtxt_node_ordering(mphtxt_t%nn(:,i), mphtxt_t%type)
          i = i+1
          if (i > mphtxt_t%nel) cycle ! Number of parameters already readed. 

        ! Local number of parameters per element
        elseif (local_nparam == -1 .and. word_count(line,' ') == 1) then
          read(line,*) local_nparam
        ! Number of parameters
        elseif (nparam == -1 .and. word_count(line,' ') == 1) then
          read(line,*) nparam
          if(nparam == 0) cycle
        ! Parameters
        elseif ((j <= nparam)) then ! .and. word_count(line,' ') == local_nparam*mphtxt_t%local_nnodes) then
          ! Skip
          j = j+1
          if (j > nparam) cycle ! Number of parameters already readed. 

        ! Number of geometric indices
        elseif (nindices == -1 .and. word_count(line,' ') == 1) then
          read(line,*) nindices
          allocate(mphtxt_t%ref(nindices))
        ! Number of parameters
        elseif (allocated(mphtxt_t%ref) .and. (k <= nindices) .and. word_count(line,' ') == 1) then
          read(line,*) mphtxt_t%ref(k)
          mphtxt_t%ref(k) = mphtxt_t%ref(k) + 1 !PMH indices starts in 1
          k = k+1
          if (k > nindices) cycle ! Number of geometric indices already readed.


        ! Number of up/down pairs
        elseif (nupdownpairs == -1 .and. word_count(line,' ') == 1) then
          read(line,*) nupdownpairs
          if(nupdownpairs == 0) exit
        ! Up/down pairs
        elseif ((nupdownpairs == 0) .or. (l <= nupdownpairs .and. word_count(line,' ') == 2)) then
          ! Skip
          l = l+1
          if (l > nupdownpairs) exit ! Number of up/down pairs already readed. Type readed.
        else
          exit
        endif
      endif
    enddo


end subroutine

function mphtxt_get_type(desc) result(res)

  character(len=*), intent(in)  :: desc
  integer                       :: res
  integer                       :: nnod, nver, lnn, lnv, lne, lnf

    nnod=0; nver=0; lnn=0; lnv=0; lne=0; lnf=0

    if(trim(desc) == 'vtx') then                     ! Node
      nnod=1; nver=1; lnn=1; lnv=1; lne=0; lnf=0
      call info('Element type: Node')
    elseif(trim(desc) == 'edg') then                 ! Edge Lagrange P1
      nnod=2; nver=2; lnn=2; lnv=2; lne=1; lnf=0
      call info('Element type: Edge lagrange P1')
    elseif(trim(desc) == 'tri') then                 ! Triangle Lagrange P1
      nnod=3; nver=3; lnn=3; lnv=3; lne=3; lnf=0
      call info('Element type: Triangle lagrange P1')
    elseif(trim(desc) == 'quad') then                ! Quadrangle Lagrange P1
      nnod=4; nver=4; lnn=4; lnv=4; lne=4; lnf=0
      call info('Element type: Quadrangle lagrange P1')
    elseif(trim(desc) == 'tet') then                 ! Tetrahedron Lagrange P1
      nnod=4; nver=4; lnn=4; lnv=4; lne=6; lnf=4
      call info('Element type: Tetrahedron lagrange P1')
    elseif(trim(desc) == 'prism') then               ! Prism Lagrange P1
      ! Prism FE not supported
      call error('Wedge lagrange P1 not supported')
    elseif(trim(desc) == 'hex') then                 ! Hexahedron Lagrange P1
      nnod=8; nver=8; lnn=8; lnv=8; lne=12; lnf=6 
      call info('Element type: Hexahedron lagrange P1')
    elseif(trim(desc) == 'edg2') then                ! Edge Lagrange P2
      nnod=3; nver=2; lnn=3; lnv=2; lne=1; lnf=0
      call info('Element type: Edge lagrange P2')
    elseif(trim(desc) == 'tri2') then                 ! Triangle Lagrange P2
      nnod=6; nver=3; lnn=6; lnv=3; lne=3; lnf=0
      call info('Element type: Triangle lagrange P2')
    elseif(trim(desc) == 'quad2') then               ! Quadrangle Lagrange P2
      nnod=9; nver=4; lnn=9; lnv=4; lne=4; lnf=0
      call info('Element type: Quadrangle lagrange P2')
    elseif(trim(desc) == 'tet2') then                ! Tetrahedron Lagrange P2
      nnod=10; nver=4; lnn=10; lnv=4; lne=6; lnf=4
      call info('Element type: Tetrahedron lagrange P2')
    elseif(trim(desc) == 'prism2') then              ! Prism Lagrange P2
      ! Quadratic prism FE not supported
      call error('Wedge lagrange P2 not supported')
    elseif(trim(desc) == 'hex2') then                ! Hexahedron Lagrange P2
      nnod=26; nver=8; lnn=16; lnv=8; lne=12; lnf=6
      call info('Element type: Hexahedron lagrange P2')
  endif

  res = check_fe(nnod==nver, lnn, lnv, lne, lnf) 

end function


subroutine mphtxt_node_ordering(el, tp)

  integer, dimension(:), intent(inout) :: el
  integer,               intent(in)    :: tp
  integer                              :: aux
  integer, dimension(:), allocatable   :: auxel

    if (tp <= 0) then
      call error('module_read_mphtxt/node_ordering # Element type not supported')
    endif

    if ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and.                       &               ! Nodes
        (FEDB(tp)%lnn == 1) .and. (FEDB(tp)%lnv == 1) .and.              &
        (FEDB(tp)%lne == 0) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in nodes

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and.                   &               ! Edge Lagrange P1
            (FEDB(tp)%lnn == 2) .and. (FEDB(tp)%lnv == 2) .and.          &
            (FEDB(tp)%lne == 1) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and.                   &               ! Triangle Lagrange P1
            (FEDB(tp)%lnn == 3) .and. (FEDB(tp)%lnv == 3) .and.          &
            (FEDB(tp)%lne == 3) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in triangles lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and.                   &               ! Quadrangle Lagrange P1
            (FEDB(tp)%lnn == 4) .and. (FEDB(tp)%lnv == 4) .and.          &
            (FEDB(tp)%lne == 4) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = MPH[1,2,4,3]
        if (size(el,1) /= FEDB(tp)%lnn) then 
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and.                   &               ! Tetrahedron Lagrange P1
            (FEDB(tp)%lnn == 4) .and. (FEDB(tp)%lnv == 4) .and.          &
            (FEDB(tp)%lne == 6) .and. (FEDB(tp)%lnf == 4)) then
        ! PMH and MPHTXT uses the same node ordering in tetrahedrons lagrange P1 

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and.                   &               ! Hexahedron Lagrange P1
            (FEDB(tp)%lnn == 8)  .and. (FEDB(tp)%lnv == 8) .and.         &
            (FEDB(tp)%lne == 12) .and. (FEDB(tp)%lnf == 6)) then
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = MPH[1,2,4,3,5,6,8,7]
        if (size(el,1) /= FEDB(tp)%lnn) then 
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(7); el(7) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and.                  &                ! Edge Lagrange P2
            (FEDB(tp)%lnn == 3) .and. (FEDB(tp)%lnv == 2) .and.          &
            (FEDB(tp)%lne == 1) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P2

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and.                  &                ! Triangle Lagrange P2
            (FEDB(tp)%lnn == 6) .and. (FEDB(tp)%lnv == 3) .and.          &
            (FEDB(tp)%lne == 3) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in triangles lagrange P2
        ! PMH[1,2,3,4,5,6] = MPH[1,2,3,4,6,5]
        if (size(el,1) /= FEDB(tp)%lnn) then 
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(6); el(6) = el(5); el(5) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and.                  &                ! Quadragle Lagrange P2
            (FEDB(tp)%lnn == 9) .and. (FEDB(tp)%lnv == 4) .and.          &
            (FEDB(tp)%lne == 4) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in quadrangles lagrange P2

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and.                  &                ! Tetrahedron Lagrange P2
            (FEDB(tp)%lnn == 10) .and. (FEDB(tp)%lnv == 4) .and.         &
            (FEDB(tp)%lne == 6)  .and. (FEDB(tp)%lnf == 4)) then
        ! PMH and MPHTXT don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10] = MPH[1,2,4,3,5,7,6,8,9,10]
        if (size(el,1) /= FEDB(tp)%lnn) then 
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(7); el(7) = el(6); el(6) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and.                   &               ! Hexahedron Lagrange P2
            (FEDB(tp)%lnn == 26) .and. (FEDB(tp)%lnv == 8) .and.          &
            (FEDB(tp)%lne == 12) .and. (FEDB(tp)%lnf == 6)) then
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26] = 
        ! MPH[1,2,4,3,5,6,8,7,9,12,13,10,14,16,22,20,23,26,27,24,11,17,15,25,19,21]

        if (size(el,1) /= 27) then 
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif

        if (allocated(auxel)) deallocate(auxel)
        allocate(auxel(size(el,1)))
        auxel(:) = el(:)

        el(1) = auxel(1); el(2) = auxel(2); el(3) = auxel(4); el(4) = auxel(3)
        el(5) = auxel(5); el(6) = auxel(6); el(7) = auxel(8); el(8) = auxel(7)
        el(9) = auxel(9); el(10) = auxel(12); el(11) = auxel(13); el(12) = auxel(10)
        el(13) = auxel(14); el(14) = auxel(16); el(15) = auxel(22); el(16) = auxel(20)
        el(17) = auxel(23); el(18) = auxel(26); el(19) = auxel(27); el(20) = auxel(24)
        el(21) = auxel(11); el(22) = auxel(17); el(23) = auxel(15); el(24) = auxel(25)
        el(25) = auxel(19); el(26) = auxel(21)

    endif


end subroutine


end module
