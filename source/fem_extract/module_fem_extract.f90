module module_fem_extract
!-----------------------------------------------------------------------
! Module for fem extractions
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Rodrigo Valina; Francisco Pena, fran.pena@usc.es
! Last update: 17/03/2011
!
! PUBLIC TYPE-BOUND PROCEDURES:
!   extract_mesh: extract a submesh from a global mesh selecting 
!   only the elements with some particular subdomain numbers
!   extract_field: extract a subfield from a global field
!   cell2node: calculate a node field form a cell one
!
! USAGE:
!   1) EXTRACT MESH: 
!      call extract_mesh(nver, mm, z, nsd, nsd0, submm, subz, globv, globel)
!
!      a) The global mesh is given by nver, mm, z, nsd
!      b) The subd. numbers to choose are given by nsd0
!      c) The submesh is returned in submm, subz
!      d) The index vectors to know the old indices of vertices and elements
!      are returned in globv, globel
!
!   2) EXTRACT FIELD:
!      call extract_field(globv, globel, v, subv, type, ncomp)
!
!      a) The index vector to know the old indices of vertices is given by globv
!      b) The index vector to know the old indices of elements is given by globel
!      c) The global field is given by v
!      d) The subfield is returned in subv
!      d) The character argument type must be 'node' or 'cell'
!      d) The optional integer argument ncomp is the component number; set to 1 by default.
!
!      Examples:
!      call extract_field(globv, globel, temp,  sub_temp,  'node')    ! scalar    node field
!      call extract_field(globv, globel, veloc, sub_veloc, 'node', 3) ! 3D vector node field
!      call extract_field(globv, globel, joule, sub_joule, 'cell')    ! scalar    cell field
!      call extract_field(globv, globel, gradV, sub_gradV, 'cell', 3) ! 3D vector cell field
!
!   2) CELL to NODE:
!      call cell2node(nver, mm, vc, vn)
!
!      a) The mesh is given by nver, mm
!      a) The cell field is given by vn
!      b) The node field is returned in vc
!-----------------------------------------------------------------------
use module_fem_extract_real
use module_fem_extract_complex
implicit none

end module
