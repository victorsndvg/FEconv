subroutine tet_mesh_order_lnn_adj_count(node_num, tet_num, tet_node, &
  adj_num, adj_row, lnn)

!*****************************************************************************80
!
!! TET_MESH_ORDER_LNN_ADJ_COUNT counts the number of nodal adjacencies.
!
!  Discussion:
!
!    Assuming that the tet mesh is to be used in a finite element
!    computation, we declare that two distinct nodes are "adjacent" if and
!    only if they are both included in some tetrahedron.
!
!    It is the purpose of this routine to determine the number of
!    such adjacency relationships.
!
!    The initial count gets only the (I,J) relationships, for which
!    node I is strictly less than node J.  This value is doubled
!    to account for symmetry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2018
!
!  Author:
!
!    John Burkardt
!    Modified by Francisco Pena
!
!  Parameters:
!
!    Input, integer(4) NODE_NUM, the number of nodes.
!
!    Input, integer(4) TET_NUM, the number of tetrahedrons.
!
!    Input, integer(4) TET_NODE(LNN,TET_NUM), the indices of
!    the nodes.
!
!    Output, integer(4) ADJ_NUM, the total number of adjacency
!    relationships,
!
!    Output, integer(4) ADJ_ROW(NODE_NUM+1), the ADJ pointer array.
!
  implicit none

  integer(4) tet_num
  integer(4) node_num
  integer(4) lnn
  integer(4) adj_num
  integer(4), dimension(node_num+1) :: adj_row
  integer(4) i
  integer(4) j
  integer(4) k
  integer(4), allocatable, dimension(:,:) :: pair
  integer(4) pair_num
  integer(4) pair_unique_num
  integer(4), dimension(lnn,tet_num) :: tet_node
!
!  Each order lnn-node tetrahedron defines lnn*(lnn-1)/2 adjacency pairs.
!
  if (.not. allocated(pair)) allocate(pair(2,lnn*(lnn-1)*tet_num/2))

  k = 0
  do i = 1, lnn-1
    do j = i + 1, lnn
      pair(1,k*tet_num+1:(k+1)*tet_num) = tet_node(i,1:tet_num)
      pair(2,k*tet_num+1:(k+1)*tet_num) = tet_node(j,1:tet_num)
      k = k + 1
    end do
  end do
!
!  Force the nodes of each pair to be listed in ascending order.
!
  pair_num = lnn*(lnn-1)*tet_num/2
  call i4col_sort2_a(2, pair_num, pair)
!
!  Rearrange the columns in ascending order.
!
  call i4col_sort_a(2, pair_num, pair)
!
!  Get the number of unique columns.
!
  call i4col_sorted_unique_count(2, pair_num, pair, pair_unique_num)
!
!  The number of adjacencies is TWICE this value, plus the number of nodes.
!
  adj_num = 2 * pair_unique_num
!
!  Now set up the ADJ_ROW counts.
!
  adj_row(1:node_num) = 0

  do k = 1, pair_num

    if(1 < k) then
      if(pair(1,k-1) == pair(1,k) .and. &
           pair(2,k-1) == pair(2,k)) then
        cycle
      end if
    end if

    i = pair(1,k)
    j = pair(2,k)

    adj_row(i) = adj_row(i) + 1
    adj_row(j) = adj_row(j) + 1

  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  return
end subroutine

subroutine tet_mesh_order_lnn_adj_set(node_num, tet_num, tet_node, &
  adj_num, adj_row, adj, lnn)

!*****************************************************************************80
!
!! TET_MESH_ORDER_LNN_ADJ_SET sets the nodal adjacency matrix.
!
!  Discussion:
!
!    A compressed format is used for the nodal adjacency matrix.
!
!    It is assumed that we know ADJ_NUM, the number of adjacency entries
!    and the ADJ_ROW array, which keeps track of the list of slots
!    in ADJ where we can store adjacency information for each row.
!
!    We essentially repeat the work of TET_MESH_ORDER_LNN_ADJ_COUNT, but
!    now we have a place to store the adjacency information.
!
!    A copy of the ADJ_ROW array is useful, as we can use it to keep track
!    of the next available entry in ADJ for adjacencies associated with
!    a given row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2018
!
!  Author:
!
!    John Burkardt
!    Modified by Francisco Pena
!
!  Parameters:
!
!    Input, integer(4) NODE_NUM, the number of nodes.
!
!    Input, integer(4) TET_NUM, the number of tetrahedrons.
!
!    Input, integer(4) TET_NODE(LNN,TET_NUM), the indices of
!    the nodes.
!
!    Input, integer(4) ADJ_NUM, the total number of adjacency
!    relationships,
!
!    Input, integer(4) ADJ_ROW(NODE_NUM+1), the ADJ pointer array.
!
!    Output, integer(4) ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer(4) adj_num
  integer(4) tet_num
  integer(4) node_num
  integer(4) lnn
  integer(4), dimension(adj_num) :: adj
  integer(4), dimension(node_num+1) :: adj_row
  integer(4), dimension(node_num+1) :: adj_row_copy
  integer(4) i
  integer(4) j
  integer(4) k
  integer(4), allocatable, dimension(:,:) :: pair
  integer(4) pair_num
  integer(4), dimension(lnn,tet_num) :: tet_node
!
!  Each order lnn-node tetrahedron defines lnn*(lnn-1)/2 adjacency pairs.
!
  if (.not. allocated(pair)) allocate(pair(2,lnn*(lnn-1)*tet_num/2))

  k = 0
  do i = 1, lnn-1
    do j = i + 1, lnn
      pair(1,k*tet_num+1:(k+1)*tet_num) = tet_node(i,1:tet_num)
      pair(2,k*tet_num+1:(k+1)*tet_num) = tet_node(j,1:tet_num)
      k = k + 1
    end do
  end do
!
!  Force the nodes of each pair to be listed in ascending order.
!
  pair_num = lnn*(lnn-1)*tet_num/2
  call i4col_sort2_a(2, pair_num, pair)
!
!  Rearrange the columns in ascending order.
!
  call i4col_sort_a(2, pair_num, pair)
!
!  Mark all entries of ADJ so we will know later if we missed one.
!
  adj(1:adj_num) = -1
!
!  Copy the ADJ_ROW array and use it to keep track of the next
!  free entry for each row.
!
  adj_row_copy(1:node_num) = adj_row(1:node_num)
!
!  Now set up the ADJ_ROW counts.
!
  do k = 1, pair_num

    if(1 < k) then
      if(pair(1,k-1) == pair(1,k) .and. &
           pair(2,k-1) == pair(2,k)) then
        cycle
      end if
    end if

    i = pair(1,k)
    j = pair(2,k)

    adj(adj_row_copy(i)) = j
    adj_row_copy(i) = adj_row_copy(i) + 1
    adj(adj_row_copy(j)) = i
    adj_row_copy(j) = adj_row_copy(j) + 1

  end do

  return
end subroutine
