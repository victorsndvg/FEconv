module module_pvd_fcnv
!-----------------------------------------------------------------------
! Module to manage PVD files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 15/05/2013
!
! PUBLIC PROCEDURES:
!   load_vtu: loads a VTU format file
!   read_vtu: reads a mesh with fields in a VTU format file
!   save_vtu_mfm: saves a mesh in a VTU format file from mfm (refs= nsd_*, nrc_*, nra_*,nrv_*)
!   save_vtu_pmh: saves a mesh and fields in a VTU format file from pmh
!   type_cell: give the associated name of FE
!-----------------------------------------------------------------------
use basicmod, only: slash, sunique
use module_vtu_fcnv
use module_pmh_fcnv
implicit none

contains

!-----------------------------------------------------------------------
! load_pvd(filename, pmh, infield): write VTU file
!-----------------------------------------------------------------------
! filename: name of a VTU file
! pmh:      PMH structure storing the piecewise mesh
! infield:  input field names
!-----------------------------------------------------------------------

subroutine load_pvd(filename, pmh, infield)
  character(len=*),             intent(in) :: filename
  type(pmh_mesh),            intent(inout) :: pmh
  character(len=*), allocatable,intent(in) :: infield(:)
  integer                                  :: iu, ios, counter, p1, p2
  logical                                  :: file_exists
  character(len=maxpath)                   :: line
  character(len=maxpath)                   :: timeatt = 'timestep='
  character(len=maxpath)                   :: fileatt = 'file='
  real(real64)                             :: timestep
  character(len=maxpath)                   :: file, path


  ! Check if file exists
  inquire(file=trim(filename), exist=file_exists)

  if(.not. file_exists) call error('Input file '//trim(filename)//' not found!')

  iu = get_unit()
  counter = 0

  open(unit=iu, file=trim(filename), form='FORMATTED', access='SEQUENTIAL', action = 'READ', iostat = ios)
  if(ios /= 0) call error("load_pvd/open #"//trim(string(ios)))

  path = filename(:index(filename,slash(), back=.true.))

  do
    read(unit=iu, fmt='(A)', iostat = ios) line
    if(ios > 0) then
      call error("load_pvd/open #"//trim(string(ios)))
    elseif(ios<0) then
      exit
    endif

    if(index(line, '<DataSet') > 0) then
      counter = counter + 1
      p1 = index(line, trim(timeatt))+len_trim(timeatt)+1
      p2 = index(line(p1:), '"')+p1-2
      read(unit=line(p1:p2), fmt=*, iostat = ios) timestep
      if(ios /= 0) call error('load_pvd/cannot convert timestep to real#'//trim(string(ios)))
      p1 = index(line, trim(fileatt))+len_trim(fileatt)+1
      p2 = index(line(p1:), '"')+p1-2
      file = trim(line(p1:p2))
      inquire(file=trim(filename), exist=file_exists)
      if(.not. file_exists) call error('PVD pointed file '//trim(filename)//' not found!')
      call info('PVD: Reading file "'//trim(path)//trim(file)//'" with timestep "'//trim(string(timestep))//'"')
      call load_vtu(trim(path)//trim(file), pmh, infield,nparam=counter,param=timestep)
    endif

  enddo

  close(unit=iu,iostat=ios)
  if(ios /= 0) call error("close_pvd/close #"//trim(string(ios)))



end subroutine

!-----------------------------------------------------------------------
! save_vtu_pmh(filename, pmh, infield, outfield, padval): write VTU file
!-----------------------------------------------------------------------
! filename: name of a VTU file
! pmh:      PMH structure storing the piecewise mesh
! infield:  input field names
! outfield: output field names
! padval:   value for padding nodes or elements without values
!-----------------------------------------------------------------------
subroutine save_pvd(filename, pmh, infield, outfield, padval)

  character(len=*),              intent(in) :: filename
  type(pmh_mesh),             intent(inout) :: pmh ! pmh_mesh
  character(len=*), allocatable, intent(in) :: infield(:)
  character(len=*), allocatable, intent(in) :: outfield(:)
  real(real64),                  intent(in) :: padval

  real(real64)                              :: p
  character(len=maxpath)                    :: fname, pvdfname, prefn, path
  integer, allocatable                      :: nshots(:),unshots(:)
  integer                                   :: i, iu, ios


  call get_num_shots(pmh, nshots)

  if(allocated(nshots)) then
    if(size(nshots,1)<1) call error("No fields found!")
    call sunique(nshots,unshots)
    if(size(unshots,1)>1) call error("Contains fields with different number of shots")
  else
    call error("No fields found!")
  endif

  prefn =  filename(index(filename,slash(), back=.true.)+1:index( filename, '.', back=.true.)-1)
  path = filename(:index(filename,slash(), back=.true.))

  call write_pvd_collection(filename, iu)
  do i=1, unshots(1)
    fname= trim(path)//trim(adjustl(prefn))//trim(string(i))//'.vtu'
    pvdfname= trim(adjustl(prefn))//trim(string(i))//'.vtu'
    call save_vtu(fname, pmh, infield, outfield, padval, nparam=i, param=p)
    call info('PVD: Saved file "'//trim(fname)//'" with timestep "'//trim(string(p))//'"')
    write(unit=iu,fmt='(A)',iostat=ios) &
      & '        <DataSet timestep="'//trim(string(p))//'" group="" part="0" file="'//trim(pvdfname)//'"/>'
    if(ios /= 0) call error("save_pvd/dataset #"//trim(string(ios)))
  enddo
  call close_pvd_collection(iu)


end subroutine

subroutine write_pvd_collection(filename, iu)
  character(len=*), intent(in) :: filename
  integer,       intent(inout) :: iu
  integer                      :: ios

  iu = get_unit()
  open(unit=iu, file=trim(filename), form='FORMATTED', access='SEQUENTIAL', action = 'WRITE', iostat = ios)
  if(ios /= 0) call error("write_pvd/open #"//trim(string(ios)))

  write(unit=iu,fmt='(A)',iostat=ios)'<?xml version="1.0"?>'
  if(ios /= 0) call error("write_pvd/version #"//trim(string(ios)))
  write(unit=iu,fmt='(A)',iostat=ios)'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
  if(ios /= 0) call error("write_pvd/vtkfile #"//trim(string(ios)))
  write(unit=iu,fmt='(A)',iostat=ios)'    <Collection>'
  if(ios /= 0) call error("write_pvd/collection #"//trim(string(ios)))

end subroutine

subroutine close_pvd_collection(iu)
  integer,       intent(inout) :: iu
  integer                      :: ios

  write(unit=iu,fmt='(A)',iostat=ios)'    </Collection>'
  if(ios /= 0) call error("close_pvd/collection #"//trim(string(ios)))
  write(unit=iu,fmt='(A)',iostat=ios)'</VTKFile>'
  if(ios /= 0) call error("close_pvd/vtkfile #"//trim(string(ios)))

  close(unit=iu,iostat=ios)
  if(ios /= 0) call error("close_pvd/close #"//trim(string(ios)))

end subroutine




end module
