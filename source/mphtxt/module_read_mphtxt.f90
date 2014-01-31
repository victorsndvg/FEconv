module module_read_mphtxt
!-----------------------------------------------------------------------
! Module for mphtxt file read
! Last update: 08/01/2014
!-----------------------------------------------------------------------
use module_COMPILER_DEPENDANT, only: real64
use module_os_dependant, only: maxpath
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
          if (version(1) /= 0 .or. version(2) /= 1) call error('mphtxt_file/header, 0.1 is the only supported version #'//string(version))
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
!    if (allocated(mphtxt_o%etypes)) deallocate(mphtxt_o%etypes)
    if (allocated(mphtxt_o%z)) deallocate(mphtxt_o%znod)
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
          if (.not. (serializable(1) == 0 .or. serializable(1) == 1)) call error('mphtxt_file/object, Only versions 0 or 1 supported #'//string(serializable))
          if (serializable(3) /= 1) call error('mphtxt_file/object, Not serializable object #'//string(serializable))
	  print*, 'serializable:', serializable
        ! Object class
        elseif (obj_class == '' .and. word_count(line,' ') == 2) then
          read(line,*) aux, obj_class
	  print*, 'obj_class:', trim(obj_class)
          if (obj_class /= 'Mesh') call error('mphtxt_file/object, Only mesh object allowed #'//trim(obj_class))
        ! Object version
        elseif (version == -1 .and. (word_count(line,' ') == 1)) then
          read(line,*) version
	  print*, 'version:', version
        ! Space dimension
        elseif ((mphtxt_o%dim == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) mphtxt_o%dim
	  print*, 'dimension:', mphtxt_o%dim
        ! Number of points
        elseif ((mphtxt_o%nnod == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) mphtxt_o%nnod
	  print*, 'nnod:', mphtxt_o%nnod
        ! Lowest mesh point index
        elseif ((offset == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) offset
	  print*, 'nnod0:', offset
          allocate(mphtxt_o%znod(mphtxt_o%dim,mphtxt_o%nnod))
        ! Point coordinates
        elseif (allocated(mphtxt_o%znod) .and. (i <= mphtxt_o%nnod) .and. (word_count(line,' ') == mphtxt_o%dim)) then
             read(line,*) (mphtxt_o%znod(k,i), k=1,mphtxt_o%dim)
	     print*, 'znod:', mphtxt_o%znod(:,i)
             i = i+1
          if (i > mphtxt_o%nnod) cycle ! All coords already readed
        elseif ((netypes == -1) .and. (word_count(line,' ') == 1)) then
          read(line,*) netypes
	  print*, 'elements:', netypes
          allocate(mphtxt_o%el(netypes))
          exit ! Number of elements already readed. Object header readed.
        else
          exit
        endif
      endif
    enddo

    ! Read object element types
    do i = 1, netypes
      print*, 'TYPE:',i
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
    print*, 'READ_MPHTXT_ETYPE'
    ! Read object header
    do 
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios /= 0) call error('mphtxt_file/object/etype, #'//trim(string(ios)))

      line = trim(line,'#')

      if (len_trim(line) /= 0) then  ! Discards empty lines
        ! FE type
        if (len_trim(fetype_name) == 0 .and. word_count(line,' ') == 2) then
          read(line,*) aux, fetype_name
          print*, 'fe type:', trim(fetype_name)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          mphtxt_t%type = get_mphtxt_eltype(fetype_name)
!          call assign_FE_type(mphtxt_t%fetype, fetype_name)
        ! Local number of nodes per element
        elseif (local_nnodes == -1 .and. word_count(line,' ') == 1) then
          read(line,*) local_nnodes
	  print*, 'local_nnodes: ', local_nnodes
        ! Number of elements
        elseif (mphtxt_t%nel == -1 .and. word_count(line,' ') == 1) then
          read(line,*) mphtxt_t%nel
	  print*, 'mphtxt_t%nel: ', mphtxt_t%nel
          allocate(mphtxt_t%nn(local_nnodes, mphtxt_t%nel))
        ! Elements
        elseif (allocated(mphtxt_t%nn) .and. (i <= mphtxt_t%nel).and. word_count(line,' ') == local_nnodes) then
          read(line,*) (mphtxt_t%nn(m,i), m=1,local_nnodes)
          mphtxt_t%nn(:,i) = mphtxt_t%nn(:,i) - offset + 1
	  if (mphtxt_t%type == QU2_P1 .or. mphtxt_t%type == QU3_P1) then
              aux = mphtxt_t%nn(4,i)
              mphtxt_t%nn(4,i) = mphtxt_t%nn(3,i)
              mphtxt_t%nn(3,i) = aux
          endif

	  if (mphtxt_t%type == ED2_P1 .or. mphtxt_t%type == ED3_P1) then
		print*, '********', mphtxt_t%nn(:,i)
          endif
	  print*, 'mphtxt_t%nn: ', mphtxt_t%nn(:,i)
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
          print*, 'nindices: ', nindices
          allocate(mphtxt_t%ref(nindices))
        ! Number of parameters
        elseif (allocated(mphtxt_t%ref) .and. (k <= nindices) .and. word_count(line,' ') == 1) then
          read(line,*) mphtxt_t%ref(k)
          mphtxt_t%ref(k) = mphtxt_t%ref(k) + 1 !PMH indices starts in 1
          print*, 'mphtxt_t%ref:', mphtxt_t%ref(k)
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


end module
