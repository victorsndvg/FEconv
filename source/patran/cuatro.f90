function cuatro(nrv1, nrv2, nrv3, nrv4, ndir) result(nrc)
!funcion cuatro: asignacion de referencias a caras (original de Inés, mas abajo)
!
! NOTA: se supone que, en la creación de los nrv se ha tenido en cuenta que:
!  - las cond. Dirichlet tienen prioridad sobre las Neumann;
!  - dentro de un tipo, las cond. con número alto, tienen prioridad sobre los bajos

implicit none
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