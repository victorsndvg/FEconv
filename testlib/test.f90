program test
use basicmod, only: real64, maxpath, error, string, lcase, file_exists, get_unit, operator(.IsNewerThan.)
use module_feconv, only: convert, pmh_mesh, pmh2mfm, save_mfm
implicit none
! MFM mesh
character(maxpath) :: meshfile
integer :: doe, dim, lnv, lne, lnf, nel, nver
integer, allocatable :: mm(:,:), nrv(:,:), nra(:,:), nrc(:,:), nsd(:)
real(real64), allocatable :: z(:,:)
! DOF
integer :: lnn, nnod
integer, allocatable :: nn(:,:)
! PMH mesh
type(pmh_mesh) :: apmh

print*, 'Simple test to check FEconv as library'
print*, '--------------------------------------'
! convert mesh to PMH using FEconv
call convert('"meshin.vtu"', outpmh=apmh)
! convert PMH to MFM
call pmh2mfm(apmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
! dimension of FE
if     (lne <= 1) then; doe = 1
elseif (lnf <= 1) then; doe = 2
else ;                  doe = 3
end if
! save mesh in MFM to disk
if (.not. file_exists('meshout.mfm') .or. ('meshin.vtu' .IsNewerThan. 'meshout.mfm')) then
  print*, 'meshout.mfm will be updated.'
  call save_mfm('meshout.mfm', get_unit(), nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
else
  print*, 'meshout.mfm will not be updated.'
end if
end program

