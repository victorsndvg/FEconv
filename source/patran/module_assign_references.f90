module module_assign_references_fcnv
!-----------------------------------------------------------------------
! Module to assign references in Patran format
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Ines Santos Atienza
!          Francisco Pena, fran.pena(at)usc.es
! Last update: 18/11/2016
!
! PUBLIC PROCEDURES:
!   load_patran: loads a Patran mesh file
!   dos: assign references to edges
!   tres: assign references to faces
!   cuatro: assign references to faces
!-----------------------------------------------------------------------
implicit none

contains

function dos(nrv1, nrv2, ndir) result(nra)
!funcion dos: asignacion de referencias a aristas (original de Inés, mas abajo)
!
! NOTA: se supone que, en la creación de los nrv se ha tenido en cuenta que:
!  - las cond. Dirichlet tienen prioridad sobre las Neumann;
!  - dentro de un tipo, las cond. con número alto, tienen prioridad sobre los bajos

!implicit none
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

function tres(nrv1, nrv2, nrv3, ndir) result(nrc)
!funcion tres: asignacion de referencias a caras (original de Inés, mas abajo)
!
! NOTA: se supone que, en la creación de los nrv se ha tenido en cuenta que:
!  - las cond. Dirichlet tienen prioridad sobre las Neumann;
!  - dentro de un tipo, las cond. con número alto, tienen prioridad sobre los bajos

!implicit none
integer, intent(in) :: nrv1, nrv2, nrv3 !numeros de referencia de los vertices de la cara
integer, intent(in) :: ndir             !numero total de condiciones Dirichlet
integer :: nrc    !numero de referencia asignado a la cara
integer :: mn, mx !minimo y maximo del los nrv

mn = min(nrv1, nrv2, nrv3); mx = max(nrv1, nrv2, nrv3)
if (mn == 0) then
  nrc = 0 !si algun nrv == 0, nrc = 0
elseif (mx <= ndir .or. mn > ndir) then
  nrc = mn !si los nrv son del mismo tipo (Dirichlet o Neumann), se toma el menor
else !hay algun vertice Neumann (pero no todos)
  nrc = nrv1 !nrv1 puede ser Neumann (o no)
  !si nrc es Dirichlet o nrv2 es Neumann y menor que nrc, se cambia a nrv2
  if (nrc <= ndir .or. (ndir < nrv2 .and. nrv2 < nrc)) nrc = nrv2
  if (nrc <= ndir .or. (ndir < nrv3 .and. nrv3 < nrc)) nrc = nrv3 !idem nrv3
endif

end function

!integer function tres(a,b,c,h)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subrutina empleada en la asignación de referencias a caras
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  implicit none
!  
!  integer a, b, c, h, min, med, max
!
!  if (a.lt.b) then
!    min=a; max=b
!  else
!    min=b; max=a
!  end if
!
!  if (c.lt.min) then
!    med=min; min=c
!  else
!	if (c.gt.max) then
!	  med=max; max=c
!	else
!	  med=c 
!	end if
!  end if
!
!
!  if (max.eq.0) then
!    tres=0
!  else
!	if (max.le.h) then
!	  tres=max
!	else
!	  if (min.gt.h) then
!	    tres=max
!	  else
!	    if (med.gt.h) then
!	      tres=min
!		else
!		  tres=med
!		end if
!	  end if
!	end if
!  end if
!
!end function

function cuatro(nrv1, nrv2, nrv3, nrv4, ndir) result(nrc)
!funcion cuatro: asignacion de referencias a caras (original de Inés, mas abajo)
!
! NOTA: se supone que, en la creación de los nrv se ha tenido en cuenta que:
!  - las cond. Dirichlet tienen prioridad sobre las Neumann;
!  - dentro de un tipo, las cond. con número alto, tienen prioridad sobre los bajos

!implicit none
integer, intent(in) :: nrv1, nrv2, nrv3, nrv4 !numeros de referencia de los vertices 
integer, intent(in) :: ndir                   !numero total de condiciones Dirichlet
integer :: nrc    !numero de referencia asignado a la cara
integer :: mn, mx !minimo y maximo del los nrv

mn = min(nrv1, nrv2, nrv3, nrv4); mx = max(nrv1, nrv2, nrv3, nrv4)
if (mn == 0) then
  nrc = 0 !si algun nrv == 0, nrc = 0
elseif (mx <= ndir .or. mn > ndir) then
  nrc = mn !si los nrv son del mismo tipo (Dirichlet o Neumann), se toma el menor
else !hay algun vertice Neumann (pero no todos)
  nrc = nrv1 !nrv1 puede ser Neumann (o no)
  !si nrc es Dirichlet o nrv2 es Neumann y menor que nrc, se cambia a nrv2
  if (nrc <= ndir .or. (ndir < nrv2 .and. nrv2 < nrc)) nrc = nrv2
  if (nrc <= ndir .or. (ndir < nrv3 .and. nrv3 < nrc)) nrc = nrv3 !idem nrv3
  if (nrc <= ndir .or. (ndir < nrv4 .and. nrv4 < nrc)) nrc = nrv4 !idem nrv4
endif

end function

!integer function cuatro(a,b,c,d,h)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Subrutina empleada en la generación de referencias a caras
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  implicit none
!  
!  integer a, b, c, d, h, min, medmin, medmax, max

!  if (a.lt.b) then
!    min=a; max=b
!  else
!    min=b; max=a
!  end if
!
!  if (c.lt.min) then
!    medmin=min; min=c
!  else
!	if (c.gt.max) then
!	  medmin=max; max=c
!	else
!	  medmin=c
!	end if
!  end if
!
!  if (d.lt.min) then
!    medmax=medmin; medmin=min; min=d
!  else
!	if (d.gt.max) then
!	  medmax=max; max=d
!	else
!	  if (d.gt.medmin) then
!	    medmax=d
!	  else
!	    medmax=medmin; medmin=d
!	  end if
!	end if
!  end if
!
!
!  if (max.eq.0) then
!    cuatro=0
!  else
!	if (max.le.h) then
!	  cuatro=max
!	else
!	  if (min.gt.h) then
!	    cuatro=max
!	  else
!	    if (medmin.gt.h) then
!	      cuatro=min
!		else
!		  if (medmax.gt.h) then
!		    cuatro=medmin
!		  else
!		    cuatro=medmax
!		  end if
!		end if
!	  end if
!	end if
!  end if
!
!end function

end module
