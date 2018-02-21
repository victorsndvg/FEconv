module module_xml_parser_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for reading XML files suitable for [OpenNum](https://sourceforge.net/projects/opennum/).  
!!
!! @note The basic rules for OpenNum-suitable XML files are:  
!! * Valid XML marks are `menu`, `submenu`, `struct` and `leaf`.  
!! * All marks must contain tag `name`.  
!! * Only mark `leaf` can store data and all its data appears inside mark `elements`.  
!! * Marks `leaf` must contain tag `totalnum`, indicating the number of stored data.  
!! See XML examples below for more information.  
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: real64
use module_os_dependant_bmod, only: maxpath
use module_alloc_bmod, only: csize, extend, alloc, dealloc
use module_convers_bmod, only: string, int, adjustlt
use module_report_bmod, only: error, info, debug
use module_files_bmod, only: get_unit
implicit none

!Constants
character(*), dimension(4), parameter, private ::  OPEN_MARK_MEMBERS = &
(/ '<menu   ', '<submenu','<struct ', '<leaf   ' /)
character(*), dimension(4), parameter, private :: CLOSE_MARK_MEMBERS = &
(/ '</menu   ', '</submenu','</struct ', '</leaf   ' /)
character(*), dimension(2), parameter, private ::  OPEN_MARK_LEAF = (/ '<leaf', '<data'/)
character(*), dimension(2), parameter, private :: CLOSE_MARK_LEAF = (/ '</leaf', '</data'/)
character(*), dimension(1), parameter, private ::  OPEN_MARK_ELEMENT = (/ '<element' /)
character(*), dimension(1), parameter, private :: CLOSE_MARK_ELEMENT = (/ '</element' /)

!Class attributes

!Private methods
private :: search_mark_once, search_mark, search_close_mark, follow_path, &
           last_part, cut_end_delimiter
private :: fread_real, fread_vreal, fread_vreal_alloc, &
           fread_complex, fread_vcomplex, fread_vcomplex_alloc, &
           fread_char, fread_vchar, fread_vchar_alloc

!Interfaces
interface fread
  !! Read data from a XML file suitable for [OpenNum](https://sourceforge.net/projects/opennum/).  
  !!
  !! @note The second argument of `fread` is a XML path; it is composed of `name` values concatenated by separators 
  !! (given as the first character of the path).  
  !!
  !! @warning In the array case, an error is raised when the argument size is smaller than the value of tag `totalnum`.  
  module procedure fread_real
  module procedure fread_vreal
  module procedure fread_complex
  module procedure fread_vcomplex
  module procedure fread_char
  module procedure fread_vchar
end interface
interface fread_alloc
  !! Read data from a XML file suitable for [OpenNum](https://sourceforge.net/projects/opennum/) using allocatable arrays.  
  !!
  !! @note The second argument of `fread_alloc` is a XML path; it is composed of `name` values concatenated by separators 
  !! (given as the first character of the path).  
  !!
  !! @warning If `realloc` is not present or it is `.true.`, the array is reallocated to save the XML data. 
  !! Otherwise, an error is raised when the argument size is smaller than the value of tag `totalnum`.  
  module procedure fread_vreal_alloc
  module procedure fread_vcomplex_alloc
  module procedure fread_vchar_alloc
end interface

contains

!-----------------------------------------------------------------------
! fopen
!-----------------------------------------------------------------------
function fopen(datxml) result(un)
!! Open a XML file.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, fopen, fread, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `real(real64) :: val`  
!! `ide = fopen()`  
!! `call fread(ide, '/Boundary conditions/Neumann Condition/Constant value', val)`  
!! `call fclose(ide)`  
!! `end program`  
!!
!! @note The XML file given in `datxml` is opened by `fopen`. If `datxml` is missing then:  
!! * If two command arguments `-xml` _filename_ are provided, then _filename_ is opened.  
!! * If the only command argument `-xml` is provided, then `local.mnu.xml` is opened.  
!! * If no command argument is provided, then `local.mnu.xml` is opened.  
character(*), intent(in), optional :: datxml !! XML file suitable for [OpenNum](https://sourceforge.net/projects/opennum/).
character(maxpath) :: xmlfile, arg1
integer :: un !associated unit number
integer :: ios, length, status

!find a valid xmlfile
xmlfile = ' '
if (present(datxml)) then
  xmlfile = datxml
else
  if (command_argument_count() == 2) then
    call get_command_argument(1, arg1, length, status)
    if (status /= 0) call error('fopen/get_command_argument, '&
     &//'the first command argument cannot be read')
    if (trim(adjustl(arg1)) == '-xml') then
      call get_command_argument(2, xmlfile, length, status)
      if (status /= 0) call error('fopen/get_command_argument, '&
       &//'the second command argument cannot be read')
    else
      call error('fopen/get_command_argument, '&
       &//'the first command argument is not recognized (must be -xml)')
    endif
  elseif (command_argument_count() == 1) then
    call get_command_argument(1, arg1, length, status)
    if (status /= 0) call error('fopen/get_command_argument, '&
     &//'the first command argument cannot be read')
    if (trim(adjustl(arg1)) == '-xml') then
      xmlfile = 'local.dat.xml'
    else
      call error('fopen/get_command_argument, '&
       &//'the first command argument is not recognized (must be -xml)')
    endif
  elseif (command_argument_count() == 0) then
    xmlfile = 'local.dat.xml'
  else
    call error('fopen/get_command_argument, '&
     &//'too many arguments')
  end if
end if
if (len_trim(xmlfile)==0) call error('fopen, filename is empty')
call info('fopen (xmlfile), '//trim(xmlfile))
! open xmlfile
un = get_unit()
open (unit=un, file=xmlfile, form='formatted', iostat=ios, &
status='old', position='rewind')
if (ios /= 0) call error('fopen, #'//trim(string(ios)))

end function

!-----------------------------------------------------------------------
! fread_real(fread)
!-----------------------------------------------------------------------
subroutine fread_real(un, path, var)
!! Read a rea64 scalar from the XML tree.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, fopen, fread, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `real(real64) :: val`  
!! `ide = fopen()`  
!! `call fread(ide, '/Boundary conditions/Neumann Condition/Constant value', val)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Constant value" totalnum="1" type="float">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`-99`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
integer,      intent(in)  :: un   !! Unit number given by `fopen`.
character(*), intent(in)  :: path !! XML path composed of `name` values.
real(real64), intent(out) :: var  !! Real64 scalar.
integer :: res, tn = 1

call follow_path(un, path, back = .true.)
if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='float', tnum=tn)) call error(trim(path)//'), not found')
!get_elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
read (un, *, iostat=res) var
if (res /= 0) call error('fread_real/read ('//trim(path)//'), #'//trim(string(res)))
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
call info('fread_real ('//trim(path)//'), '//trim(string(var)))
end subroutine

!-----------------------------------------------------------------------
! fread_vreal (fread)
!-----------------------------------------------------------------------
subroutine fread_vreal(un, path, var)
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, fopen, fread, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `real(real64) :: arr(5)`  
!! `ide = fopen()`  
!! `call fread(ide, '/Boundary conditions/Neumann Condition/Interval', arr)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Interval" totalnum="2" type="float">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`0.01 0.05`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
integer,      intent(in)    :: un     !! Unit number given by `fopen`.
character(*), intent(in)    :: path   !! XML path composed of `name` values.
real(real64), intent(inout) :: var(:) !! Real64 array.
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='float', tnum=tn)) call error(trim(path)//'), not found')
if (tn > csize(dbl1=var, d=1)) call error('fread_vreal ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than array size '//trim(string(csize(dbl1=var, d=1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vreal/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  call info('fread_vreal ('//trim(path)//'), '//trim(string(var(i))))
enddo
end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_complex(un, path, var)

integer,          intent(in)  :: un
character(len=*), intent(in)  :: path
character(len=128)            :: tempstring
complex(real64),  intent(out) :: var
integer :: res, tn = 1

call follow_path(un, path, back = .true.)
if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='complex', tnum=tn)) call error(trim(path)//'), not found')
!get_elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
read (un, *, iostat=res) var
if (res /= 0) call error('fread_complex/read ('//trim(path)//'), #'//trim(string(res)))
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
write(tempstring,*) var
call info('fread_complex ('//trim(path)//'), '//trim(tempstring))

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vcomplex(un, path, var)

integer,          intent(in)  :: un
character(len=*),           intent(in)  :: path
character(len=128)                      :: tempstring
complex(real64), dimension(:), intent(out) :: var
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='complex', tnum=tn)) call error(trim(path)//'), not found')
if (tn > size(var,1)) call error('fread_vcomplex ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vcomplex/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  write(tempstring,*) var(i)
  call info('fread_vcomplex ('//trim(path)//'), '//trim(tempstring))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread: read data from xml file
!-----------------------------------------------------------------------
subroutine fread_vcomplex_alloc(un, path, var, realloc)

integer,                    intent(in)  :: un
character(len=*),           intent(in)  :: path
character(len=128)                      :: tempstring
complex(real64), dimension(:), intent(inout), allocatable :: var
logical,                    intent(in), optional :: realloc
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='complex', tnum=tn)) call error(trim(path)//'), not found')
if (present(realloc)) then; if (realloc) then
  if (allocated(var)) deallocate(var)
endif; endif
if (.not. allocated(var)) allocate(var(tn))
if (tn > size(var,1)) call error('fread_vcomplex_alloc ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than expected '//trim(string(size(var,1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vcomplex_alloc/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  write(tempstring,*) var(i)
  call info('fread_vcomplex_alloc ('//trim(path)//'), '//trim(tempstring))
enddo

end subroutine

!-----------------------------------------------------------------------
! fread_char (fread)
!-----------------------------------------------------------------------
subroutine fread_char(un, path, var)
!! Read a character scalar from the XML tree.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: fopen, fread, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `character(20) :: val`  
!! `ide = fopen()`  
!! `call fread(ide, '/Boundary conditions/Neumann Condition/Label', val)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Label" totalnum="1" type="charlist">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Upper side`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
integer,      intent(in)  :: un   !! Unit number given by `fopen`.
character(*), intent(in)  :: path !! XML path composed of `name` values.
character(*), intent(out) :: var  !! Character scalar.
integer :: res, tn = 1

call follow_path(un, path, back = .true.)
if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), tnum=tn)) & !typ can be diverse
  call error(trim(path)//'), not found')
!get_elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
read (un, '(a)', iostat=res) var
var = adjustlt(var)
if (res /= 0) call error('fread_char/read ('//trim(path)//'), #'//trim(string(res)))
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
call info('fread_char ('//trim(path)//'), '//trim(var))
end subroutine

!-----------------------------------------------------------------------
! fread_vchar (fread)
!-----------------------------------------------------------------------
subroutine fread_vchar(un, path, var)
!! __Example:__  
!! `program test`  
!! `use basicmod, only: fopen, fread, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `character(20) :: arr(5)`  
!! `ide = fopen()`  
!! `call fread(ide, '/Boundary conditions/Neumann Condition/Surfaces', arr)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Surfaces" totalnum="2" type="charlist">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Plate 1`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Bar 2`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
integer,      intent(in)    :: un     !! Unit number given by `fopen`.
character(*), intent(in)    :: path   !! XML path composed of `name` values.
character(*), intent(inout) :: var(:) !! Character array.
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not. search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), tnum=tn)) & !typ can be diverse
  call error('fread_vchar ('//trim(path)//'), not found')
if (tn > csize(char1=var, d=1)) call error('fread_vchar ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than array size '//trim(string(csize(char1=var, d=1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, '(a'//trim(string(len(var(1))))//')', iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vreal/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  var(i) = adjustlt(var(i))
  call info('fread_vchar ('//trim(path)//'), '//trim(var(i)))
enddo
end subroutine

!-----------------------------------------------------------------------
! fread_vreal_alloc (fread_alloc)
!-----------------------------------------------------------------------
subroutine fread_vreal_alloc(un, path, var, realloc)
!! Read numeric data from a XML file using allocatable arrays.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, fopen, fread_alloc, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `real(real64), allocatable :: arr(:)`  
!! `ide = fopen()`  
!! `call fread_alloc(ide, '/Boundary conditions/Neumann Condition/Interval', arr)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Interval" totalnum="2" type="float">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`0.01 0.05`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
integer,      intent(in)    :: un     !! Unit number given by `fopen`.
character(*), intent(in)    :: path   !! XML path composed of `name` values.
real(real64), intent(inout), allocatable :: var(:)  !! Real64 allocatable array.
logical,      intent(in), optional       :: realloc !! Whether to realloc `var` if necessary.
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not.search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), &
typ='float', tnum=tn)) call error(trim(path)//'), not found')
if (present(realloc)) then; if (realloc) call alloc(var, tn); endif
if (.not. present(realloc)) call alloc(var, tn) ! by default, realloc = .true.
if (tn > csize(dbl1=var, d=1)) call error('fread_vreal_alloc ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than array size '//trim(string(csize(dbl1=var, d=1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, *, iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vreal/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  call info('fread_vreal_alloc ('//trim(path)//'), '//trim(string(var(i))))
enddo
end subroutine

!-----------------------------------------------------------------------
! fread_vchar_alloc (fread_alloc)
!-----------------------------------------------------------------------
subroutine fread_vchar_alloc(un, path, var, realloc)
!! Read alphanumeric data from a XML file using allocatable arrays.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: fopen, fread_alloc, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `character(20), allocatable :: arr(:)`  
!! `ide = fopen()`  
!! `call fread_alloc(ide, '/Boundary conditions/Neumann Condition/Surfaces', arr)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Surfaces" totalnum="2" type="charlist">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Plate 1`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Bar 2`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
integer,      intent(in)    :: un     !! Unit number given by `fopen`.
character(*), intent(in)    :: path   !! XML path composed of `name` values.
character(*), intent(inout), allocatable :: var(:)  !! Character allocatable array.
logical,      intent(in), optional       :: realloc !! Whether to realloc `var` if necessary.
integer :: res, tn, i

call follow_path(un, path, back = .true.)
!get the total number
tn = -1; if (.not. search_mark_once(un, path, OPEN_MARK_LEAF, name = last_part(path), tnum=tn)) & !typ can be diverse
  call error('fread_vchar_alloc ('//trim(path)//'), not found')
if (present(realloc)) then; if (realloc) call alloc(var, tn); endif
if (.not. present(realloc)) call alloc(var, tn) ! by default, realloc = .true.
if (tn > csize(char1=var, d=1)) call error('fread_vreal_alloc ('//trim(path)//'), found totalnum '//&
trim(string(tn))//' is bigger than array size '//trim(string(csize(char1=var, d=1))))
!get elements
if (.not.search_mark_once(un, path, OPEN_MARK_ELEMENT)) call error(trim(path)//', not found')
if (tn>0) then
    read (un, '(a'//trim(string(len(var(1))))//')', iostat=res) var(1:tn)
    if (res /= 0) call error('fread_vchar/read ('//trim(path)//'), #'//trim(string(res)))
endif
if (.not.search_mark_once(un, path,CLOSE_MARK_ELEMENT)) call error(trim(path)//', not found')
do i = 1, tn
  var(i) = adjustlt(var(i))
  call info('fread_vchar_alloc ('//trim(path)//'), '//trim(var(i)))
enddo
end subroutine

!-----------------------------------------------------------------------
! flist
!-----------------------------------------------------------------------
subroutine flist(un, path, var)
!! Read `name` values from the children of a given XML path.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: fopen, flist, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `character(20), allocatable :: lst(:)`  
!! `ide = fopen()`  
!! `call flist(ide, '/Boundary conditions/Neumann Condition/', lst)`  
!! `call fclose(ide)`  
!! `end program`  
!! A valid XML input for the previous program would be the following _local.dat.xml_ file:  
!! `<?xml version='1.0' encoding='iso-8859-15'?>`  
!! `<data>`  
!! &nbsp;`<menu name="Boundary conditions">`  
!! &nbsp;&nbsp;`<submenu name="Neumann Condition">`  
!! &nbsp;&nbsp;&nbsp;`<leaf name="Surfaces" totalnum="2" type="charlist">`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`<elements>`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Plate 1`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`Bar 2`  
!! &nbsp;&nbsp;&nbsp;&nbsp;`</elements>`  
!! &nbsp;&nbsp;&nbsp;`</leaf>`  
!! &nbsp;&nbsp;`</submenu>`  
!! &nbsp;`</menu>`  
!! `</data>`  
!! @warning The array is always reallocated to save all the values.  
integer,      intent(in)    :: un     !! Unit number given by `fopen`.
character(*), intent(in)    :: path   !! XML path composed of `name` values.
character(*), intent(inout), allocatable :: var(:)  !! Character allocatable array.
integer :: p, n, res
logical :: close_mark_found
character(10*maxpath), allocatable :: tmp(:)
character(10*maxpath) :: line

n = 0
!initial allocation of tmp
if (.not. allocated(tmp)) call alloc(tmp, 10)
!follow the path (advancing the last line)
call follow_path(un, path)
!check whether a close mark is reached
close_mark_found = search_mark_once(un, path, CLOSE_MARK_MEMBERS, back = .true.)
!loop to catch members
do while (.not. close_mark_found)
  !searchs a new member
  call search_mark(un, path, OPEN_MARK_MEMBERS, back=.true.)
  !stores the name
  read (un, fmt='(A)', iostat=res) line
  if (res /= 0) call error('flist/read ('//trim(path)//'), #'//trim(string(res)))
  p = index(line, 'name=')
  n = n + 1; if (n > csize(char1=tmp, d=1)) call extend(tmp, n, fit=.false.)
  tmp(n) = trim(string(cut_end_delimiter(line(p+5:),'>')))
  !searchs the end of the member
  call search_close_mark(un, path)
  !check whether a close mark is reached
  close_mark_found = search_mark_once(un, path, CLOSE_MARK_MEMBERS, back = .true.)
enddo
!advance one line in file
close_mark_found = search_mark_once(un, path, CLOSE_MARK_MEMBERS)
!final result
if (allocated(var)) call dealloc(var)
call alloc(var, n)
var(1:n) = tmp(1:n)
end subroutine

!-----------------------------------------------------------------------
! fclose
!-----------------------------------------------------------------------
subroutine fclose(un)
!! Close a XML file.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: real64, fopen, fread, fclose`  
!! `implicit none`  
!! `integer :: ide`  
!! `real(real64) :: val`  
!! `ide = fopen()`  
!! `call fread(ide, '/Boundary conditions/Neumann Condition/Constant value', val)`  
!! `call fclose(ide)`  
!! `end program`  
integer, intent(in) :: un !! Unit number given by `fopen`.
integer :: ios
close(unit=un, iostat=ios)
if (ios /= 0) call error('read_xml/close, #'//trim(string(ios)))
end subroutine

!-----------------------------------------------------------------------
! PRIVATE  PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! search_mark_once: searchs a mark only once
! RETURN: .true. if the mark is found
!         .false. otherwise
!-----------------------------------------------------------------------
recursive function search_mark_once(un, path, marks, name, typ, tnum, advance, back) result(res)
integer,          intent(in) :: un
character(*), intent(in) :: path
character(*), intent(in), dimension(:) :: marks
character(*), intent(in), optional :: name, typ
integer,          intent(inout), optional :: tnum
logical, intent(in), optional :: back, advance
logical :: res
integer :: ios, i, p
character(len=10*MAXPATH) :: line

res = .false.
!read a line
read (un, fmt='(A)', iostat=ios) line
if (ios /= 0) call error('search_mark_once/read ('//trim(path)//'), #'//trim(string(ios)))
!backspace
if (present(back)) then
  if (back) then
    backspace(unit=un, iostat=ios)
    if (ios /= 0) call error('search_mark_once/backspace ('//trim(path)//'), #'//trim(string(ios)))
  endif
endif
do i = 1, size(marks, 1)
  if (index(line, trim(marks(i))) > 0) then
    !check name
    if (present(name)) then
      p = index(line, 'name=')
      if (trim(name) /= trim(string(cut_end_delimiter(line(p+5:),'>')))) then
        if (present(advance)) then
          if (advance) call search_close_mark(un, path)
        endif
        cycle
      endif
    endif
    !check type
    if (present(typ)) then
      p = index(line, 'type=')
      if (trim(typ) /= trim(string(cut_end_delimiter(line(p+5:),'>')))) then
        if (present(advance)) then
          if (advance) call search_close_mark(un, path)
        endif
        cycle
      endif
    endif
    !check totalnum
    if (present(tnum)) then
      p = index(line, 'totalnum=')
      if (tnum > 0) then
        if (tnum /= int(string(cut_end_delimiter(line(p+9:),'>')))) then
          if (present(advance)) then
            if (advance) call search_close_mark(un, path)
          endif
          cycle
        endif
       else
        tnum = int(string(cut_end_delimiter(line(p+9:),'>')))
      endif
    endif
    !the mark, name, type and/or tnum matches
    res = .true.; return
  endif
enddo

end function

!-----------------------------------------------------------------------
! search_mark: searchs a mark
!-----------------------------------------------------------------------
subroutine search_mark(un, path, marks, name, typ, tnum, back, advance)
integer,          intent(in)  :: un
character(*), intent(in) :: path
character(*), intent(in), dimension(:) :: marks
character(*), intent(in), optional :: name, typ
integer,          intent(inout), optional :: tnum
logical, intent(in), optional :: back, advance
integer :: ios

do
  if (search_mark_once(un, path, marks, name, typ, tnum, advance, back=.false.)) then
    !backspace
    if (present(back)) then; if (back) then
      backspace(unit=un, iostat=ios)
      if (ios /= 0) call error('search_mark/backspace, ('//trim(path)//') #'//trim(string(ios)))
    endif; endif
    return
  endif
  !ends the loop if a close mark is found
  if (search_mark_once(un, path, CLOSE_MARK_MEMBERS, back = .true.)) exit
enddo
!mark not found
call error('search_mark ('//trim(path)//'), not found')

end subroutine

!-----------------------------------------------------------------------
! search_close_mark: search a close mark
!-----------------------------------------------------------------------
recursive subroutine search_close_mark(un, path)

integer,           intent(in)  :: un
character(*),  intent(in) :: path
integer :: n

n = 1 !number of open marks
do while (n > 0)
  if (search_mark_once(un, path, OPEN_MARK_MEMBERS, back=.true.)) n = n + 1
  if (search_mark_once(un, path, CLOSE_MARK_MEMBERS)) n = n - 1
enddo

end subroutine

!-----------------------------------------------------------------------
! follow_path: follow the path
!-----------------------------------------------------------------------
subroutine follow_path(un, path, back)
integer,           intent(in) :: un
character(*),      intent(in) :: path
logical, optional, intent(in) :: back
character(len=len(path)) :: lpath, parte
character(len=1) :: separador
integer :: p, res

rewind(un)
separador = path(1:1)
lpath = path(2:)
parte = lpath
do while (len_trim(parte) > 0)
  p = index(lpath, separador)
  if (p > 0) then
    parte = lpath(:p-1); lpath = lpath(p+1:)
  else
    parte = lpath;       lpath = ' '
  endif
  call search_mark(un, path, OPEN_MARK_MEMBERS, parte, advance=.true.)
  parte = lpath
enddo
if (present(back)) then; if (back) then
  backspace(unit=un, iostat=res)
  if (res /= 0) call error('follow_path/backspace ('//trim(path)//'), #'//trim(string(res)))
endif; endif
end subroutine

!-----------------------------------------------------------------------
! last_part: extracts the last part of a path
!-----------------------------------------------------------------------
function last_part(path) result(res)

character(*),  intent(in) :: path
character(len=len(path)) :: res
character(len=1) :: separador
integer :: p

separador = path(1:1)
p = index(path, separador, back=.true.)
res = path(p+1:)

end function

!-----------------------------------------------------------------------
! cut_end_delimiter: cuts an end delimiter
!-----------------------------------------------------------------------
function cut_end_delimiter(str, delimiter) result(res)
character(*), intent(in) :: str, delimiter
character(len=len(str)) :: res
integer :: p

res = str
p = index(str, delimiter)
if  (p > 0) res = str(:p-1)

end function

end module
