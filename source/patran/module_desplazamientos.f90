module module_desplazamientos_fcnv
!-----------------------------------------------------------------------
! Modulo para guardar las condiciones sobre el desplazamiento
! Last update: 27/07/2009
! Programmer: fran.pena@usc.es
!
! ATRIBUTOS PRIVADOS DE CLASE:
!   dof: matrix 6 x n de grados de libertad
!   val: matrix 6 x n de valores del desplazamiento
!   (cada fila de 'dof' y 'val' indica los grados y valores de la condicion n-esima)
!   nod: para cada n, nodos asociados a la condicion n-esima
!
! METODOS PUBLICOS:
!   set_SPC: para cada linea con una condicion SPC, almacena (dof, val, nod)
!   assign_SPC: asigna los numeros de referencia Dirichlet a cada nodo
!-----------------------------------------------------------------------
 use basicmod, only: string
 use module_ALLOC_int_alloc_r2_fcnv
 use module_ALLOC_real_r2_fcnv
 use module_ALLOC_log_r2_fcnv
 implicit none

!Variables
integer, private                              :: ncond = 0  !numero total de condiciones
logical, private, dimension(:,:), allocatable :: dof        !grados de libertad
real,    private, dimension(:,:), allocatable :: val        !valores del desplazamiento
type(int_alloc_r2), private                   :: nod        !nodos asociados a la condicion
integer, private                              :: nucond = 0 !numero total de condiciones unicas
integer, private, dimension(:),   allocatable :: ucond      !indices (para nod) de las condiciones unicas

contains

!-----------------------------------------------------------------------
! set_SPC: para cada linea con una condicion SPC, almacena (dof, val, nod)
!
! Lectura del cammpo SPC: MD Nastran 2006. Quick Reference Guide (p.2467)
! SPC: (A8)
! SID: Identification number of the single-point constraint set. (I8)
! G1:  Grid or scalar point identification number. (I8)
! C1:  Component number. (0 < Integer < 6; up to six Unique Integers)
! D1: Value of enforced motion for all degrees-of-freedom (Real, 8 pos.)
!-----------------------------------------------------------------------
subroutine set_SPC(iu)

  integer, intent(in) :: iu !identificador del fichero
  integer :: SID, G1
  character(len=8) :: SPC, C1, D1
  logical, dimension(6) :: newdof
  real,    dimension(6) :: newval, tmp
  integer :: res, i

  newdof = .false.
  newval = 0.
! lectura del registro
  read(unit = iu, fmt = '(A8,I8,I8,A8,A8)', iostat = res) SPC, SID, G1, C1, D1
  if (res /= 0) call info('(module_desplazamientos/set_SPC) Unable to read record')
! construcci�n de newdof, newval
  do i = 1, 6
    if (index(C1, trim(string(i))) > 0) newdof(i) = .true. !'i' aparece en C1
  end do
  where (newdof) newval = real(D1)
  if (ncond > 0) then !ya existen condiciones previas en nod
    if (any(nod%row(ncond)%col == G1)) then !se encontro G1 en la ultima condicion
      where (newdof); tmp = 1.; elsewhere; tmp = 0.; end where
      if (.not. any(dof(:,ncond) .and. (/(.false., i=1,maxloc(tmp,1)-1), &
                                         (.true.,  i=maxloc(tmp,1),6)/))) then
        !la ultima condicion, en la cual aparece el nodo, es compatible con newdof:
        !(ninguno de los dof(:,ncond) a partir del primer .true. en newdof esta activado)
        call set(dof, newdof, col=ncond, fit_row=.true., add=.true.) !se anhaden los nuevos dof
        where (newdof); tmp = newval; elsewhere; tmp = 0.; end where !in F95, pack...
        call set(val, tmp, col=ncond, fit_row=.true., add=.true.) !se anhaden los nuevos val
        return !salimos sin ejecutar el resto del codigo
      end if
    end if
  end if
  !o no se encontro G1 o son incompatibles
  ncond = ncond + 1
  call set(dof, newdof, col=ncond, fit_row=.true.)
  call set(val, newval, col=ncond, fit_row=.true.)
  call set(nod, G1, row=ncond, fit_col=.true.)

end subroutine

!-----------------------------------------------------------------------
!   assign_SPC: asigna los numeros de referencia Dirichlet a cada nodo
!-----------------------------------------------------------------------
subroutine assign_SPC(nnrv, nD)

  integer, dimension(*) :: nnrv !numeros de referencia (por nodo)
  integer, intent(in)   :: nD   !numero total de condiciones Dirichlet (ya asignadas)
  integer :: i, j, ni, nj, k

! eliminacion de condiciones duplicadas
  do i = ncond, 2, -1 !bucle en condiciones, de la ultima a la primera
    do j = 1, i - 1 !bucle en condiciones previas a la i-esima
      if (all(dof(:,i).eqv.dof(:,j)) .and. maxval(abs(val(:,i)-val(:,j)))<epsilon(1.)) then
!       dof y val coinciden, se agrupan los nodos en j
        ni = size(nod%row(i)%col, 1)
        nj = size(nod%row(j)%col, 1)
        do k = 1, ni !adicion de las condiciones i-�simas a las condiciones j-�simas
          call set(nod, nod%row(i)%col(k), row=j, col=nj+k) !por rapidez, no se ajusta en cols
        end do
        call dealloc(nod%row(i)%col) !eliminaci�n de las condiciones i-�simas
        call reduce(nod%row(j)%col, ni + nj) !al final se reduce en cols
!        ncond = ncond - 1 !eliminacion de la ultima condicion
        exit
      end if
    end do
  end do

! detecci�n de condiciones unificadas (se construye nucond, ucond)
  do i = 1, ncond
    if (allocated(nod%row(i)%col)) then !guardamos las condiciones que tienen nodos asociados
      call set(ucond, i, nucond + 1) !no es necesario ajustar
      nucond = nucond + 1
    end if
  end do

! asignacion de condiciones
  do i = 1, nucond
    do j = 1, size(nod%row(ucond(i))%col, 1)
      nnrv(nod%row(ucond(i))%col(j)) = nD + i !de esta forma, se asigna el mayor valor
    end do
  end do

end subroutine

end module
