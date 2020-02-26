module module_read_mphtxt_fcnv

!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Collaborator: Luis Perez
! Last update: 26/02/2020
!
! PUBLIC PROCEDURES:
! read_mphtxt_header: read the header of the MPHTXT file
! read_mphtxt_object: read a pieze of the MPHTXT mesh
! read_mphtxt_etype:  read a element group of the MPHTXT mesh
!-----------------------------------------------------------------------

use basicmod
use module_pmh_fcnv
use module_utils_mphtxt_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! read_mphtxt_header(iu, mphtxt_m): read mphtxt file header
!-----------------------------------------------------------------------
! iu:       unit number of the MPHTXT file
! mphtxt_m: PMH structure for store the mphtxt mesh
!-----------------------------------------------------------------------

subroutine read_mphtxt_header(iu, mphtxt_m)

  integer, intent(in) :: iu ! Unit number for mphtxtfile
  type(pmh_mesh), intent(inout) :: mphtxt_m
  integer,dimension(2) :: version ! Version of mphtxtfile
  character(len=MAXPATH) :: line
  integer :: aux, i, j, ios
  integer :: ntags, ntypes
  character(len=MAXPATH),dimension(:),allocatable :: tags, types ! Tags and typess

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
      if (ios /= 0) call error('read_mphtxt/header, #'//trim(string(ios)))

      line = trim(line,'#')

      if (len_trim(line) /= 0) then ! Discards empty lines
        ! File version
        if (version(1) == -1 .and. version(2) == -1 .and. word_count(line,' ') == 2) then
          read(line,*) version(1),version(2)
          if (version(1) /= 0 .or. version(2) /= 1) call error('read_mphtxt/header, 0.1 is the only supported version #'//&
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
! read_mphtxt_object(iu, mphtxt_o): read a object from mphtxt file
!-----------------------------------------------------------------------
! iu:       unit number of the MPHTXT file
! mphtxt_o: PMH pieze to store the objects contined in the MPHTXT file
!-----------------------------------------------------------------------

subroutine read_mphtxt_object(iu, mphtxt_o)

  integer, intent(in) :: iu ! Unit number for mphtxtfile
  type(piece), intent(inout) :: mphtxt_o ! mphtxt mesh
  integer, dimension(3) :: serializable ! Serializable option
  character(len=MAXPATH) :: obj_class ! Object class. Must by 'Mesh'.
  integer :: version ! Object version
  character(len=MAXPATH) :: line
  integer :: aux, i, k, ios
  integer :: netypes ! Number of element types
  integer :: offset ! lowest number of node

    ! Variables initialization
    if (allocated(mphtxt_o%z)) deallocate(mphtxt_o%z)
! if (allocated(mphtxt_o%etypes)) deallocate(mphtxt_o%etypes)
    mphtxt_o%dim = -1
    mphtxt_o%nnod = -1
    offset = -1
    netypes = -1
    serializable = (/-1,-1,-1/)
    obj_class = ''
    version = -1
    i = 1

    ! Read object header
    do
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios /= 0) call error('mphtxt_file/object, #'//trim(string(ios)))

      line = trim(line,'#')

      if (len_trim(line) /= 0) then ! Discards empty lines
        ! Serializable object
        if (serializable(1) == -1 .and. serializable(2) == -1 .and. serializable(3) == -1 .and. word_count(line,' ') == 3) then
          read(line,*) serializable(1), serializable(2), serializable(3)
          if (.not. (serializable(1) == 0 .or. serializable(1) == 1)) &
            &call error('mphtxt_file/object, Only object versions 0 or 1 supported #'//string(serializable))
          if (serializable(3) /= 1) call error('mphtxt_file/object, Not serializable object #'//string(serializable))
        ! Object class
        elseif (obj_class == '' .and. word_count(line,' ') == 2) then
          read(line,*) aux, obj_class
          if (obj_class /= 'Mesh') call error('mphtxt_file/object, Only mesh object allowed #'//trim(obj_class))
        ! Object version
        elseif (version == -1 .and. (word_count(line,' ') == 1)) then
          read(line,*) version
           if (version /= 2 .and. version /= 4) &
            &call error('mphtxt_file/object, Only Mesh versions 2 and 4 are supported #'//string(serializable))
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
      call read_mphtxt_etype(iu, mphtxt_o%el(i), offset, version)
    enddo

end subroutine


!-----------------------------------------------------------------------
! read_mphtxt_etype(iu, mphtxt_t, offset): read a element type from mphtxt file
!-----------------------------------------------------------------------
! iu:       unit number of the MPHTXT file
! mphtxt_t: PMH element group to store the element types contined in the MPHTXT file
! offset:   lowest number of nodes
!-----------------------------------------------------------------------

subroutine read_mphtxt_etype(iu, mphtxt_t, offset, version)

  integer, intent(in) :: iu ! Unit number for mphtxtfile
  type(elgroup), intent(inout) :: mphtxt_t ! mphtxt mesh
  integer, intent(in) :: offset ! Lowest number of node
  integer, intent(in) :: version !Mesh version
  character(len=MAXPATH) :: fetype_name ! Object class. Must by 'Mesh'.
  character(len=MAXPATH) :: line
  integer :: aux, i, j,k ,l, m, ios
  integer :: local_nparam, nparam, local_nnodes, nupdownpairs, nindices

    ! Variables initialization
    if (allocated(mphtxt_t%nn)) deallocate(mphtxt_t%nn)
    if (allocated(mphtxt_t%ref)) deallocate(mphtxt_t%ref)
    fetype_name = ''
    mphtxt_t%nel = -1
    local_nnodes = -1
    nindices = -1
    local_nparam = -1
    nparam = -1
    nupdownpairs = -1
    i = 1
    j = 1
    k = 1
    l = 1

    ! Read object header
    do
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios == iostat_end) then ! EOF found
        if (len_trim(fetype_name) > 0 .and. local_nnodes /= -1 .and. mphtxt_t%nel /= -1 .and. allocated(mphtxt_t%nn)) then
          ! Basic element information was read before EOF
          if (.not. allocated(mphtxt_t%ref)) then
            ! Reference information was not read; set references to 0 and exit
            allocate(mphtxt_t%ref(mphtxt_t%nel))
            mphtxt_t%ref = 0
          end if
          return
        else
          call error('(module_read_mphtxt/read_mphtxt_etype) Basic element information is missing; EOF found')
        end if
      elseif (ios /= 0) then ! An error different from EOF was found
        call error('mphtxt_file/object/etype, #'//trim(string(ios)))
      end if
      line = trim(line,'#')
      if (len_trim(line) /= 0) then ! Discards empty lines
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
          call pmh_node_ordering(mphtxt_t%nn(:,i), mphtxt_t%type)
          i = i+1
          if (i > mphtxt_t%nel) cycle ! Number of parameters already readed.

        ! Local number of parameters per element
        elseif (local_nparam == -1 .and. word_count(line,' ') == 1 .and. version == 2) then
          read(line,*) local_nparam
        ! Number of parameters
        elseif (nparam == -1 .and. word_count(line,' ') == 1 .and. version == 2) then
          read(line,*) nparam
          if(nparam == 0) cycle
        ! Parameters
        elseif (j <= nparam .and. version == 2) then ! .and. word_count(line,' ') == local_nparam*mphtxt_t%local_nnodes) then
          ! Skip
          j = j+1
          if (j > nparam) cycle ! Number of parameters already readed.

        ! Number of geometric indices, that is, references
        elseif (nindices == -1 .and. word_count(line,' ') == 1) then
          read(line,*) nindices
          allocate(mphtxt_t%ref(nindices))
        ! Number of parameters
        elseif (allocated(mphtxt_t%ref) .and. (k <= nindices) .and. word_count(line,' ') == 1) then
          read(line,*) mphtxt_t%ref(k)
          k = k+1
          if (k > nindices) then
            if (version == 2) then
              cycle ! Number of geometric indices already readed.
            elseif (version == 4) then
              exit
            else
               call error('(module_read_mphtxt/read_mphtxt_etype) Only mesh version 2 and 4 are supported.') 
            end if
          end if  
        ! Number of up/down pairs
        elseif (nupdownpairs == -1 .and. word_count(line,' ') == 1 .and. version == 2) then
          read(line,*) nupdownpairs
          if(nupdownpairs == 0) exit
        ! Up/down pairs
        elseif ((nupdownpairs == 0) .or. (l <= nupdownpairs .and. word_count(line,' ') == 2) .and. version == 2) then
          ! Skip
          l = l+1
          if (l > nupdownpairs) exit ! Number of up/down pairs already readed. Type readed.
        elseif (version == 2) then
          exit
        endif
      endif
    enddo

end subroutine



end module
