module module_fuerzas_fcnv
!-----------------------------------------------------------------------
! Modulo para guardar las condiciones sobre la fuerza
! Last update: 30/07/2009
! Programmer: fran.pena@usc.es
!
! ATRIBUTOS PRIVADOS DE CLASE:
!   val: matrix 3 x n de valores de la fuerza
!   (cada fila de 'val' indica los valores de la condicion n-esima)
!   nod: para cada n, nodos asociados a la condicion n-esima
!
! METODOS PUBLICOS:
!   set_FORCE: para cada linea con una condicion SPC, almacena (val, nod)
!   assign_FORCE: asigna los numeros de referencia Neumann a cada nodo
!-----------------------------------------------------------------------
 use basicmod, only: string
 use module_ALLOC_int_alloc_r2_fcnv
 use module_ALLOC_real_r2_fcnv
 implicit none

!Variables
integer, private                              :: ncond = 0  !numero total de condiciones
real,    private, dimension(:,:), allocatable :: val        !valores de la fuerza
type(int_alloc_r2), private                   :: nod        !nodos asociados a la condicion

contains

!-----------------------------------------------------------------------
! set_SPC: para cada linea con una entrada FORCE, almacena (val, nod)
!
! Lectura del cammpo SPC: MD Nastran 2006. Quick Reference Guide (p.1550)
! FORCE: (A8)
! SID: Load set identification number. (I8)
! G: Grid point identification number. (I8)
! CID: Coordinate system identification number. (I8)
! F: Scale factor. (Real, 8 pos.)
! N1, N2, N3: Components of a vector [...] (Real, 8 pos.)
!-----------------------------------------------------------------------
subroutine set_FORCE(iu)

  integer, intent(in) :: iu !identificador del fichero
  integer :: SID, G, CID
  character(len=8) :: FORCE, F, N1, N2, N3
  real, dimension(3) :: newval
  integer :: res, i

  newval = 0.
! lectura del registro
  read(unit=iu, fmt='(A8,I8,I8,I8,A8,A8,A8,A8)', iostat=res) FORCE, SID, G, CID, F, N1, N2, N3
  if (res /= 0) call info('(module_desplazamientos/set_FORCE) Unable to read record')
  newval = (/ real(N1), real(N2), real(N3) /) !construcci�n de newval

  do i = 1, ncond
    if (maxval(abs(val(:,i)-newval)) < epsilon(1.)) then !condicion ya guardada (en i)
      call set(nod, G, row=i, fit_col=.true.)
      exit
    end if
  end do
  if (i > ncond) then !la condicion no estaba guardada
    ncond = ncond + 1
    call set(val, newval, col=ncond, fit_row=.true.)
    call set(nod, G, row=ncond, fit_col=.true.)
  end if

end subroutine

!-----------------------------------------------------------------------
!   assign_FORCE: asigna los numeros de referencia Neumann a cada nodo
!-----------------------------------------------------------------------
subroutine assign_FORCE(nnrv, nD)

  integer, dimension(:) :: nnrv !numeros de referencia (por nodo)
  integer, intent(in)   :: nD   !numero total de condiciones Dirichlet (ya asignadas)
  integer :: i, j, n

  do i = 1, ncond
    do j = 1, size(nod%row(i)%col, 1)
      n = nod%row(i)%col(j)
!     si no tiene Dirichlet, se le asigna el mayor valor de Neumann
      if (nnrv(n) <=0 .or. nD < nnrv(n)) nnrv(n) = nD + i
    end do
  end do

end subroutine

end module
