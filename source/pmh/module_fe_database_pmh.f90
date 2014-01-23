module module_fe_database_pmh
!-----------------------------------------------------------------------
! Module to manage the finite elements database for PMH format
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 19/01/2014
!
! PUBLIC PROCEDURES:
!   check_fe_param: given the finite element parameters, return the index of its position in the database
!   check_fe_type: given the finite element description, return the index of its position in the database
!
! REMARKS:
!   To check the type of a mesh written in another format:
!      it = check_fe_param(...)
!      if (it == 0) call error('...: element type not implemented for PMH format')
!      pmesh%pc(...)%el(...)%type = fedb%type(it)
!      ...
!
!   To check whether a PMH mesh can be written in another format:
!      use module_alloc
!      (The integer array valid_fe is created using check_fe_param)
!      it = check_fe_type(pmesh%pc(...)%el(...)%type)
!      if (find_first(valid_fe, it) == 0) call error(...: element type ... cannot be written')
!-----------------------------------------------------------------------
use module_os_dependant, only: maxpath
use module_report, only: error

implicit none

!Constants
!Element types
integer, parameter :: NODE2 = 1, NODE3 = 2, &                                   !nodes
ED2_P1 =  3, ED3_P1 =  4, ED2_P2 =  5, ED3_P2 =  6, &                           !edges 
TR2_P1 =  7, TR3_P1 =  8, TR2_P2 =  8, TR3_P2 = 10, TR2_RT = 11, TR3_RT = 12, & !triangles
QU2_P1 = 13, QU3_P1 = 14, QU2_P2 = 15, QU3_P2 = 16, QU2_RT = 17, QU3_RT = 18, & !quadrangles
TET_P1 = 19, TET_P2 = 20, TET_RT = 21, TET_ND = 22, &                           !tetrahedra
HEX_P1 = 23, HEX_P2 = 24, HEX_RT = 25, HEX_ND = 26, &                           !hexahedra
WED_P1 = 27, WED_P2 = 28, WED_RT = 29, WED_ND = 30                              !wedges

!Types
type :: fe_db_pmh
  character(len=46) :: type = ' '             !description
  logical           :: nver_eq_nnod = .false. !is the total number of vertices equal to the total number of nodes?
  integer           :: dim  = 0               !space dimension
  integer           :: tdim = 0               !topological dimension (0 for nodes, 1 for edges, etc.)
  integer           :: lnn  = 0               !local number of nodes
  integer           :: lnv  = 0               !local number of vertices
  integer           :: lne  = 0               !local number of edges
  integer           :: lnf  = 0               !local number of faces
end type

!Constants                                     char(46)                d td   n  v   e  f
type(fe_db_pmh), parameter :: FEDB(30) = [ &
  fe_db_pmh('Node, dimension 2                             ', .true.,  2, 0,  1, 1,  0, 0), &
  fe_db_pmh('Node, dimension 3                             ', .true.,  3, 0,  1, 1,  0, 0), &
  fe_db_pmh('Edge, Lagrange P1, dimension 2                ', .true.,  2, 1,  2, 2,  1, 0), &
  fe_db_pmh('Edge, Lagrange P1, dimension 3                ', .true.,  3, 1,  2, 2,  1, 0), &
  fe_db_pmh('Edge, Lagrange P2, dimension 2                ', .false., 2, 1,  3, 2,  1, 0), &
  fe_db_pmh('Edge, Lagrange P2, dimension 3                ', .false., 3, 1,  3, 2,  1, 0), &
  fe_db_pmh('Triangle, Lagrange P1, dimension 2            ', .true.,  2, 2,  3, 3,  3, 1), &
  fe_db_pmh('Triangle, Lagrange P1, dimension 3            ', .true.,  3, 2,  3, 3,  3, 1), &
  fe_db_pmh('Triangle, Lagrange P2, dimension 2            ', .false., 2, 2,  6, 3,  3, 1), &
  fe_db_pmh('Triangle, Lagrange P2, dimension 3            ', .false., 3, 2,  6, 3,  3, 1), &
  fe_db_pmh('Triangle, Raviart-Thomas (edge), dimension 2  ', .false., 2, 2,  3, 3,  3, 1), &
  fe_db_pmh('Triangle, Raviart-Thomas (edge), dimension 3  ', .false., 3, 2,  3, 3,  3, 1), &
  fe_db_pmh('Quadrangle, Lagrange P1, dimension 2          ', .true.,  2, 2,  4, 4,  4, 1), &
  fe_db_pmh('Quadrangle, Lagrange P1, dimension 3          ', .true.,  3, 2,  4, 4,  4, 1), &
  fe_db_pmh('Quadrangle, Lagrange P2, dimension 2          ', .false., 2, 2,  8, 4,  4, 1), &
  fe_db_pmh('Quadrangle, Lagrange P2, dimension 3          ', .false., 3, 2,  8, 4,  4, 1), &
  fe_db_pmh('Quadrangle, Raviart-Thomas (edge), dimension 2', .false., 2, 2,  4, 4,  4, 1), &
  fe_db_pmh('Quadrangle, Raviart-Thomas (edge), dimension 3', .false., 3, 2,  4, 4,  4, 1), &
  fe_db_pmh('Tetrahedron, Lagrange P1                      ', .true.,  3, 3,  4, 4,  6, 4), &
  fe_db_pmh('Tetrahedron, Lagrange P2                      ', .false., 3, 3, 10, 4,  6, 4), &
  fe_db_pmh('Tetrahedron, Raviart-Thomas (face)            ', .false., 3, 3,  4, 4,  6, 4), &
  fe_db_pmh('Tetrahedron, Nedelec (edge)                   ', .false., 3, 3,  6, 4,  6, 4), &
  fe_db_pmh('Hexahedron, Lagrange P1                       ', .true.,  3, 3,  8, 8, 12, 6), &
  fe_db_pmh('Hexahedron, Lagrange P2                       ', .false., 3, 3, 16, 8, 12, 6), &
  fe_db_pmh('Hexahedron, Raviart-Thomas (face)             ', .false., 3, 3,  6, 8, 12, 6), &
  fe_db_pmh('Hexahedron, Nedelec (edge)                    ', .false., 3, 3, 12, 8, 12, 6), &
  fe_db_pmh('Wedge, Lagrange P1                            ', .true.,  3, 3,  6, 6,  9, 5), &
  fe_db_pmh('Wedge, Lagrange P2                            ', .false., 3, 3, 12, 6,  9, 5), &
  fe_db_pmh('Wedge, Raviart-Thomas (face)                  ', .false., 3, 3,  5, 6,  9, 5), &
  fe_db_pmh('Wedge, Nedelec (edge)                         ', .false., 3, 3,  9, 6,  9, 5)]


contains

!-----------------------------------------------------------------------
! check_fe_param: given the finite element parameters, return the index of its position in the database
!-----------------------------------------------------------------------
function check_fe_param(nnod, nver, dim, lnn, lnv, lne, lnf) result(res)
  integer, intent(in) :: nnod, nver, dim, lnn, lnv, lne, lnf
  integer :: res, i

  res = 0
  do i = 1, size(FEDB,1)
    if (((nver==nnod).eqv.FEDB(i)%nver_eq_nnod) .and. dim==FEDB(i)%dim .and. lnn==FEDB(i)%lnn .and. &
                            lnv==FEDB(i)%lnv .and. lne==FEDB(i)%lne .and. lnf==FEDB(i)%lnf) then
      res = i
      return
    end if
  enddo
end function

!-----------------------------------------------------------------------
!   check_fe_type: given the finite element description, return the index of its position in the database
!-----------------------------------------------------------------------
function check_fe_type(desc) result(res)
character(*), intent(in) :: desc
integer :: res, i

res = 0
do i = 1, size(FEDB,1)
  if (desc == FEDB(i)%type) then
    res = i
    return
  end if
enddo
end function


function get_fe_data(n) result (res)
  integer, intent(in) :: n
  type(fe_db_pmh)     :: res

  if (n<=size(FEDB,1)) then
    res = FEDB(n)
  else
    call error("get_fe_data, Wrong finite element type")
  endif
end function


function get_mphtxt_eltype(descriptor) result(res)
  character(len=MAXPATH),              intent(in)    :: descriptor  ! FE descriptor
  integer                                           :: res         ! Element type

  if (descriptor == 'vtx') then
    res = NODE2
  elseif (descriptor == 'edg') then
    res = ED2_P1
  elseif (descriptor == 'tri') then
    res = TR2_P1
  elseif (descriptor == 'tet') then
    res = TET_P1
  elseif (descriptor == 'quad') then
    res = QU2_P1
  elseif (descriptor == 'prism') then
    res = WED_P1
  elseif (descriptor == 'hex') then
    res = HEX_P1
  elseif (descriptor == 'edg2') then
    res = ED2_P2
  elseif (descriptor == 'tr2') then
    res = TR2_P2
  elseif (descriptor == 'tet2') then
    res = TET_P2
  elseif (descriptor == 'quad2') then
    res = QU2_P2
  elseif (descriptor == 'prism2') then
    res = WED_P2
  elseif (descriptor == 'hex2') then
    res = HEX_P2
  else
    call error('mphtxt/get_mphtxt_eltype, Finite element type not allowed #'//trim(descriptor))
  endif
  


end function


end module
