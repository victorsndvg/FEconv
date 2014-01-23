function dos(nrv1, nrv2, ndir) result(nra)
!funcion dos: asignacion de referencias a aristas (original de Inés, mas abajo)
!
! NOTA: se supone que, en la creación de los nrv se ha tenido en cuenta que:
!  - las cond. Dirichlet tienen prioridad sobre las Neumann;
!  - dentro de un tipo, las cond. con número alto, tienen prioridad sobre los bajos

implicit none
integer, intent(in) :: nrv1, nrv2 !numeros de referencia de los vertices de la arista
integer, intent(in) :: ndir       !numero total de condiciones Dirichlet
integer :: nra    !numero de referencia asignado a la arista
integer :: mn, mx !minimo y maximo del los nrv

mn = min(nrv1, nrv2); mx = max(nrv1, nrv2)
if (mn == 0) then
  nra = 0 !si algun nrv == 0, nra = 0
elseif (mx <= ndir .or. mn > ndir) then
  nra = mn !si los nrv son del mismo tipo(Dirichlet o Neumann), se toma el menor
elseif (nrv1 > ndir) then
  nra = nrv1 !nrv1 es Neumann
else
  nra = nrv2 !nrv2 es Neumann
endif

end function

!CODIGO ORIGINAL de INES:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subrutina empleado en la asignación de referencias a aristas
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  implicit none
!  
!  integer a, b, h, min, max
!
!  if (a.lt.b) then
!    min=a; max=b
!  else
!    min=b; max=a
!  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(min.eq.0)then
!  dos=max
!else
!  dos=min
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  if (min.eq.0) then
!    dos=0
!  else
!	if (max.le.h) then
!	  dos=min
!	else
!	  if (min.ge.h) then
!	    dos=min
!	  else
!	    dos=max
!	  end if
!	end if
!  end if

!end function
