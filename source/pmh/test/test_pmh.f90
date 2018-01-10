program test_pmh
!-----------------------------------------------------------------------
! Utility to convert between several mesh and FE field formats
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: see variable 'last_update'
!-----------------------------------------------------------------------
use module_compiler_dependant, only:real64
use module_pmh_fcnv
use module_mfm_fcnv
use module_vtu_fcnv
use module_unv_fcnv
implicit none

type(pmh_mesh) :: pmh, pmh2, pmh3

!Variables for MFM format
integer                                   :: nel  = 0 !global number of elements
integer                                   :: nnod = 0 !global number of nodes
integer                                   :: nver = 0 !global number of vertices
integer                                   :: dim  = 0 !space dimension
integer                                   :: lnn  = 0 !local number of nodes
integer                                   :: lnv  = 0 !local number of vertices
integer                                   :: lne  = 0 !local number of edges
integer                                   :: lnf  = 0 !local number of faces
integer,      allocatable, dimension(:,:) :: nn       !nodes index array
integer,      allocatable, dimension(:,:) :: mm       !vertices index array
integer,      allocatable, dimension(:,:) :: nrv      !vertices reference array
integer,      allocatable, dimension(:,:) :: nra      !edge reference array
integer,      allocatable, dimension(:,:) :: nrc      !face reference array
real(real64), allocatable, dimension(:,:) :: z        !vertices coordinates array
integer,      allocatable, dimension(:)   :: nsd      !subdomain index array


allocate(pmh%pc(2))
!Pieza 1
pmh%pc(1)%nnod = 9
pmh%pc(1)%dim  = 3
allocate(pmh%pc(1)%z(3,9))
pmh%pc(1)%z(1:3, 1) = real([0,  0,  0], real64)
pmh%pc(1)%z(1:3, 2) = real([1,  0,  0], real64)
pmh%pc(1)%z(1:3, 3) = real([1,  1,  0], real64)
pmh%pc(1)%z(1:3, 4) = real([0,  1,  0], real64)
pmh%pc(1)%z(1:3, 5) = real([0.5,0., 0.], real64)
pmh%pc(1)%z(1:3, 6) = real([1., 0.5,0.], real64)
pmh%pc(1)%z(1:3, 7) = real([0.5,1., 0.], real64)
pmh%pc(1)%z(1:3, 8) = real([0., 0.5,0.], real64)
pmh%pc(1)%z(1:3, 9) = real([0.5,0.5,0.], real64)
allocate(pmh%pc(1)%el(2))
!Grupo 1
pmh%pc(1)%el(1)%nel  = 2
pmh%pc(1)%el(1)%type = 5 !tria P2
allocate(pmh%pc(1)%el(1)%nn(6,2))
pmh%pc(1)%el(1)%nn(1:6,1) = [1,2,3,5,6,9]
pmh%pc(1)%el(1)%nn(1:6,2) = [1,3,4,9,7,8]
allocate(pmh%pc(1)%el(1)%ref(2))
pmh%pc(1)%el(1)%ref = [1,2]
!Grupo 2
pmh%pc(1)%el(2)%nel  = 1
pmh%pc(1)%el(2)%type = 3 !edge P2
allocate(pmh%pc(1)%el(2)%nn(3,1))
pmh%pc(1)%el(2)%nn(1:3,1) = [1,4,8]
allocate(pmh%pc(1)%el(2)%ref(1))
pmh%pc(1)%el(2)%ref = [11]

!Pieza 2
pmh%pc(2)%nnod = 9
pmh%pc(2)%dim  = 3
allocate(pmh%pc(2)%z(3,9))
pmh%pc(2)%z(1:3, 1) = real([1,  0,   0], real64)
pmh%pc(2)%z(1:3, 2) = real([2,  0,   0], real64)
pmh%pc(2)%z(1:3, 3) = real([2,  1,   0], real64)
pmh%pc(2)%z(1:3, 4) = real([1,  1,   0], real64)
pmh%pc(2)%z(1:3, 5) = real([1.5,0.,  0.], real64)
pmh%pc(2)%z(1:3, 6) = real([2.,  0.5,0.], real64)
pmh%pc(2)%z(1:3, 7) = real([1.5,1.,  0.], real64)
pmh%pc(2)%z(1:3, 8) = real([1., 0.5, 0.], real64)
pmh%pc(2)%z(1:3, 9) = real([1.5,0.5, 0.], real64)
!Grupo 1
allocate(pmh%pc(2)%el(2))
pmh%pc(2)%el(1)%nel  = 2
pmh%pc(2)%el(1)%type = 5 !tria P2
allocate(pmh%pc(2)%el(1)%nn(6,2))
pmh%pc(2)%el(1)%nn(1:6,1) = [1,2,3,5,6,9]
pmh%pc(2)%el(1)%nn(1:6,2) = [1,3,4,9,7,8]
allocate(pmh%pc(2)%el(1)%ref(2))
pmh%pc(2)%el(1)%ref = [1,2]
!Grupo 2
pmh%pc(2)%el(2)%nel  = 1
pmh%pc(2)%el(2)%type = 3 !edge P2
allocate(pmh%pc(2)%el(2)%nn(3,1))
pmh%pc(2)%el(2)%nn(1:3,1) = [2,3,6]
allocate(pmh%pc(2)%el(2)%ref(1))
pmh%pc(2)%el(2)%ref = [12]

call build_vertices(pmh)
call save_pmh('test.pmh', 10, pmh)

call save_unv('test.unv', 10, pmh)

!print*,'test: ', allocated(pmh%pc(1)%z), allocated(pmh%pc(1)%z)
call pmh2mfm(pmh, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)

!print*, 'nel', nel
!print*, 'nnod', nnod
!print*, 'nver', nver
!print*, 'dim', dim
!print*, 'lnn', lnn
!print*, 'lnv', lnv
!print*, 'lne', lne
!print*, 'lnf', lnf
!print*, 'nn', nn
!print*, 'mm', mm
!print*, 'nrc', nrc
!print*, 'nra', nra
!print*, 'nrv', nrv
!print*, 'z', z
!print*, 'nsd', nsd

!print*,'a'
!volver a pasarlo a PMH: pmh2
call mfm2pmh(nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd, pmh2)

!volver a pasarlo a MFM
nel = 0; nnod = 0; nver = 0; dim = 0; lnn = 0; lnv = 0; lne = 0; lnf = 0
call pmh2mfm(pmh2, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)

call save_mfm('prueba.mfm', 10, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
!print*,'b'
call save_vtu('prueba.vtu', nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
!print*,'c'


!Cuadrangulos
allocate(pmh3%pc(1))
!Pieza 1
pmh3%pc(1)%nnod = 13
pmh3%pc(1)%dim  =  3
allocate(pmh3%pc(1)%z(3,13))
pmh3%pc(1)%z(1:3, 1) = real([0,  0,  0], real64)
pmh3%pc(1)%z(1:3, 2) = real([1,  0,  0], real64)
pmh3%pc(1)%z(1:3, 3) = real([1,  1,  0], real64)
pmh3%pc(1)%z(1:3, 4) = real([0,  1,  0], real64)
pmh3%pc(1)%z(1:3, 5) = real([0.5,0., 0.], real64)
pmh3%pc(1)%z(1:3, 6) = real([1., 0.5,0.], real64)
pmh3%pc(1)%z(1:3, 7) = real([0.5,1., 0.], real64)
pmh3%pc(1)%z(1:3, 8) = real([0., 0.5,0.], real64)
pmh3%pc(1)%z(1:3, 9) = real([2,  0,  0], real64)
pmh3%pc(1)%z(1:3,10) = real([2,  1,  0], real64)
pmh3%pc(1)%z(1:3,11) = real([1.5,0. ,0.], real64)
pmh3%pc(1)%z(1:3,12) = real([2., 0.5,0.], real64)
pmh3%pc(1)%z(1:3,13) = real([0.5,1., 0.], real64)
allocate(pmh3%pc(1)%el(2))
!Grupo 1
pmh3%pc(1)%el(1)%nel  = 2
pmh3%pc(1)%el(1)%type = 7 !quad P1
allocate(pmh3%pc(1)%el(1)%nn(4,2))
pmh3%pc(1)%el(1)%nn(1:4,1) = [1,2, 4,3] !bad orientation! (is corrected by the program...)
pmh3%pc(1)%el(1)%nn(1:4,2) = [2,9,10,3]
allocate(pmh3%pc(1)%el(1)%ref(2))
pmh3%pc(1)%el(1)%ref = [1,2]
!Grupo 2
pmh3%pc(1)%el(2)%nel  = 2
pmh3%pc(1)%el(2)%type = 2 !edge P1
allocate(pmh3%pc(1)%el(2)%nn(3,2))
pmh3%pc(1)%el(2)%nn(1:2,1) = [1, 4]
pmh3%pc(1)%el(2)%nn(1:2,2) = [9,10]
allocate(pmh3%pc(1)%el(2)%ref(2))
pmh3%pc(1)%el(2)%ref = [11, 12]

call build_vertices(pmh3)
deallocate(nn, mm, nrv, nra, z, nsd)
nver = 0; nel = 0; nnod = 0
call pmh2mfm(pmh3, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
call save_vtu('prueba-quad.vtu', nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)

end program
