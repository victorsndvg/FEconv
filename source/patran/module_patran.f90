module module_patran_fcnv
!-----------------------------------------------------------------------
! Module to manage Patran mesh files
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Ines Santos Atienza
!          Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 31/05/2013
!
! PUBLIC PROCEDURES:
!   load_patran: loads a Patran mesh file
!
! REMARKS:
!   The file must have "MD Nastran input file" with extension .bdf.
!   The entries must be saved "Small Field Format" (10 fields of 8 characters)
!   References are constructed from entries SPC, FORCE of the "Bulk section"
!   For more information, consult "MD Nastran 2006 Quick Reference Guide"
!-----------------------------------------------------------------------
use basicmod, only: real64, iostat_end, det
use module_desplazamientos_fcnv
use module_fuerzas_fcnv
!use module_caras_interiores
!use module_RECONVXX
use module_groups_fcnv
use module_assign_references_fcnv, only: dos, tres, cuatro
implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                     !!
  !!   INTERFACE DE CONVERSION  PATRAN-MALLA MODULEF     !!
  !!                                                     !!
  !!                                                     !!
  !!                  Patran 2003 r2                     !!
  !!                                                     !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!  In�s Santos Atienza         25-10-2004

  !!!!!!!!!!  Departamento de Matem�tica Aplicada. USC


  !!!!!!!!!!  Modificaciones
  !
  !    27-10-2004        zn,z y zz, definidas como     real,allocatable
  !                      pasan a definirse como     doble precision,allocatable
  !

  !!!!!!!!!!  Modificaciones (Fran Pena) Julio de 2009
  !
  ! Las modificaciones se refieren a las entradas SPC, FORCE de la "Bulk section"
  ! de un fichero con formato "MD Nastran input file". Se supone que las entradas
  ! han sido almacenadas en "Small Field Format" (10 campos de 8 caracteres)
  ! Para m�s informacion, ver "MD Nastran 2006 Quick Reference Guide"
  !
  ! (el numero de linea es aproximado, los nuevos cambios lo modifican)
  !
  ! Linea 49: nombres de fichero, de character(20) a character(255)
  !
  ! Linea 103: cambio en la apertura del fichero, apertura mas general
  !
  ! Linea 404: rewind(3), permite que la lectura de coordenadas y elementos
  ! pueda hacerse en cualquier orden
  !
  ! Lineas 517 - 553: inclusion de palabra clave SPC para condiciones Dirichlet
  ! La estructura es: SPC     |SID     |G1      |C1      |D1
  !                           |        |(mm)    |(123456)|valor cond. Dirichlet
  !
  ! Linea 570: En el caso CTRIA3 y CTETRA, si el determinante de la matriz de la
  ! transformaci�n de paso al elemento de referencia es negativo, se intercambian
  ! los dos ultimos nodos
  !
  ! Lineas 560 - 622: Las entradas SPC se gestionan con module_desplazamientos.
  ! Esto permite distinguir las condiciones por los grados de libertad fijados
  ! (de 1 a 6) y los valores de desplazamiento asignados. En una entrada SPC
  ! se pueden fijar todos los grados de libertad de un nodo, o bien repartirlos
  ! en varias entradas (por ejemplo, en una primera entrada fijar los 3 primeros y,
  ! en la siguiente entrada, los 3 ultimos).
  !
  ! Supondremos que cuando se fijan los grados de libertad de una condicion y nodo
  ! concretos, se hace de forma CONTIGUA en el fichero y ORDENADA por grados de
  ! libertad (para fijar los grados 1 y 3, primero se fija el 1 y luego el 3).
  ! Este sistema es el que usa el conversor de Hypermesh a Nastran input files.
  ! De este modo, no conocemos el numero de condiciones de contorno Dirichlet hasta que
  ! se han leido todas. Por eso primero llamamos a set_SPC, que almacena todas las
  ! condiciones de entradas SPC y luego llamamos a assign:SPC, que elimina condiciones
  ! duplicadas y las asigna a nnrv
  !
  ! Ademas, el n. de ahora el n. de ref. guardado en Dirichlet no es el ultimo, sino
  ! el m�s grande. Esto nos permite en los procedimientos de asignacion de nra y nrc
  ! priorizar, para cada tipo de condicion (Neumann o Dirichlet) el numero m�s bajo.
  !
  ! Las entradas SPC1 se siguen leyendo con el codigo de Ines, pero ahora se les asigna
  ! la condici�n Dirichlet mayor.
  !
  ! Linea 644: los valores de SPC se vuelcan en nnrv
  !
  ! Linea 661: ahora el n. de ref. guardado en Neumann no es el ultimo, sino el m�s grande
  !
  ! Linea 685 - 706: Las entradas FORCE se gestionan con module_fuerzas. Esto permite
  ! distinguir las condiciones por la fuerza aplicada, no por su SID.
  ! Ademas, el n. de ref. guardado en Neumann no es el ultimo, sino el m�s grande
  !
  !!!!!!!!!!  Modificaciones (Fran Pena) Abril de 2010
  !
  ! Se refieren a la construccion de nrc en el caso de los tetraedros
  !
  !!!!!!!!!!  Modificaciones (Fran Pena) Noviembre de 2011
  !
  ! Linea 819: Modificacion para la lectura de nra mediante CBEAM, para CTRIA3 y CTRIA6

!interface
!  function dos(nrv1, nrv2, ndir) result(nra)
!  implicit none
!  integer, intent(in) :: nrv1, nrv2
!  integer, intent(in) :: ndir
!  integer :: nra
!  end function
!end interface
!interface
!  function tres(nrv1, nrv2, nrv3, ndir) result(nrc)
!  implicit none
!  integer, intent(in) :: nrv1, nrv2, nrv3
!  integer, intent(in) :: ndir
!  integer :: nrc
!  end function
!end interface
!interface
!  function cuatro(nrv1, nrv2, nrv3, nrv4, ndir) result(nrc)
!  implicit none
!  integer, intent(in) :: nrv1, nrv2, nrv3, nrv4
!  integer, intent(in) :: ndir
!  integer :: nrc
!  end function
!end interface

private :: order, swap

contains

!-----------------------------------------------------------------------
! load_patran: loads a Patran mesh file
!-----------------------------------------------------------------------
subroutine load_patran(filein, iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
character(*),              intent(in)  :: filein
integer,                   intent(in)  :: iu
integer,                   intent(out) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
integer,      allocatable, intent(out) :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64), allocatable, intent(out) :: z(:,:)
integer :: ios
character(4) :: grid
character(8) :: elto = ' ', topo, gridd, spc, through, force
!character(len=1) :: caras_interiores !Modificacion_Fran
!character(255) :: filein, fileout !Modificacion_Fran
integer :: i, j, k, kk, h, g, itopo, nmat, nodo(8), status, hd, npe, &
vpe, ape, cpe, nnnod, io !form
integer, allocatable :: m(:,:), mmver(:), mmnod(:), nnrv(:)
real(real64), allocatable :: zn(:,:), zz(:,:)
integer :: tmp
character(maxpath) :: str !Modificacion_Fran

!Modificacion_Fran
type topology
  character(8) :: name  !topology name
  integer      :: itopo !topology index (as used in this module)
  integer      :: fedim !FE dimension (<= space dimension)
  integer      :: lnn   !local number of nodes
  integer      :: lnv   !local number of vertices
  integer      :: lne   !local number of edges
  integer      :: lnf   !local number of faces
end type
type(topology) :: t(6) = [topology('CTRIA3  ', 1, 2, 3, 3,  3, 1), &
                          topology('CTRIA6  ', 2, 2, 6, 3,  3, 1), &
                          topology('CQUAD4  ', 4, 2, 4, 4,  4, 1), &
                          topology('CQUAD8  ', 5, 2, 8, 4,  4, 1), &
                          topology('CTETRA  ', 7, 3, 4, 4,  6, 4), &
                          topology('CHEXA   ', 9, 3, 8, 8, 12, 6)]
character(len=8) :: typef
integer :: ind, ref, pos, fedim = 0
integer, dimension(3) :: vnum
integer, dimension(3,4) :: face  = reshape([1,3,2, 1,4,3, 1,2,4, 2,3,4],[3,4])
integer, dimension(2,3) :: edge  = reshape([1,2, 2,3, 3,1],[2,3])
integer, dimension(2,6) :: edge2 = reshape([1,4, 4,2, 2,5, 5,3, 3,6, 6,1],[2,6])
!Fin_Modificacion_Fran

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! grid, elto, topo, gridd, spc, through, force: variables de caracteres que leen las
  !        correspondientes palabras del archivo bdf
  ! filein:   nombre del archivo de entrada (hay que escribir el '.bdf')
  ! fileout:  nombre del archivo de salida

  ! i,j,k,kk,h,g: variables �ndice

  ! itopo: integer identificador de la topolog�a de los elementos
  ! form: indica si el archivo de salida es con o sin formato

  ! nel: n�mero de elementos
  ! nnod: n�mero de nodos Patran
  ! nnnod: n�mero de nodos Modulef
  ! nver: n�mero de v�rtices
  ! nmat: n�mero de materiales distintos

  ! npe: nodos por elemento
  ! vpe: v�rtices por elemento
  ! ape: aristas por elemento
  ! cpe: caras por elemento

  ! m: matriz auxiliar para leer conectividades de v�rtices
  ! mm: matriz de conectividades de v�rtices
  ! nn: matriz de conectividades de los nodos
  ! mmver: correspondencia entre los nodos Patran y los v�rtices Modulef
  ! mmnod: correspondencia entre los nodos Patran y los nodos Modulef
  ! nnrv: vector auxiliar de referencias de v�rtices (numeraci�n Patran)
  ! nrv: referencia de los v�rtices
  ! nra: referencia de las aristas
  ! nrc: referencia de las caras
  ! nsd: n�mero de subdominio

  ! zn: matriz auxiliar para leer las coordenadas de los v�rtices
  ! z: matriz de coordenadas de los v�rtices
  ! zz: matriz de coordenadas de los nodos

  ! hd: n�mero de condiciones Dirichlet
  ! nodo: variable auxiliar para leer las referencias de los v�rtices
  ! status: variable auxiliar que se emplea en la lectura de las ref. de los v�rtices

  ! dos: tres, cuatro: output de las subrutinas con el mismo nombre

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  write(*,*)'Nombre del archivo .bdf'
!  read(*,*)filein !con .bdf!
  !Modificacion_Fran
  !open(3,file=filein,form='formatted')
  !rewind(3)
  !elto='ninguno'
  !read(3,100,iostat=io)elto
!  open(unit=3, file=filein, form='formatted', status='old', position='rewind', iostat=io)
!  if (io == 0) then
!  end if
  !Fin_Modificacion_Fran
!  do while(io.ne.0)
!    close(3)
!    write(*,*)'No se encuentra el archivo especificado. Intentalo de nuevo.'
!    write(*,*)' Recuerda escribir el nombre y la extension del archivo.'
!   read(*,*)filein
!   open(3,file=filein,form='formatted')
!      rewind(3)
!      elto='ninguno'
!      read(3,100,iostat=io)
!  end do

!  write(*,*)'Numero de nodos'
!  read(*,*)nnod !Se puede obtener de Patran. Nodos consecutivos.
!  do while(nnod<1)
!    write(*,*)' Numero de nodos no valido. Intentalo de nuevo'
!    read(*,*)nnod
!  end do
!  write(*,*)'Numero de elementos'
!  read(*,*)nel     !Se puede obtener de Patran. Elementos consecutivos.
!  do while(nel<1)
!    write(*,*)' Numero de elementos no valido. Intentalo de nuevo'
!    read(*,*)nel
!  end do
!  write(*,*)'Selecciona la topologia de elementos empleada'
!  write(*,*)'  1. Tria3        4. Quad4'
!  write(*,*)'  2. Tria6        5. Quad8'
!  write(*,*)'  3. Tria6*       6. Quad8*'
!  write(*,*)
!  write(*,*)'  7. Tet4         8. Wedge6         9. Hex8'
!  read(*,*)itopo
!  do while((itopo.lt.1).or.(itopo.gt.9))
!    write(*,*)' Opcion no valida. Intentalo de nuevo.'
!    read(*,*)itopo
!  end do
!  write(*,*)'Numero de distintos materiales'
!  read(*,*)nmat
!  do while(nmat.lt.1)
!    write(*,*)' Numero de distintos materiales no valido. Intentalo de nuevo'
!    read(*,*)nmat
!  end do
!  write(*,*)'Nombre del archivo de salida'
!  read(*,*)fileout
!  write(*,*)'Elige una opcion para el fichero de salida'
!  write(*,*)'  1. con formato'
!  write(*,*)'  2. sin formato'
!  read(*,*)form
!  do while((form.lt.1).or.(form.gt.2))
!    write(*,*)' Opcion no valida. Intentalo de nuevo'
!    read(*,*)form
!  end do
  !Modificacion_Fran
!  write(*,*) 'La malla puede tener referencias interiores para las caras (s/n)'
!  read(*,*) caras_interiores
!Fin_Modificacion_Fran
!  rewind(iu)

!open file
open (unit=iu, file=filein, form='formatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('(module_patran/load_patran) unable to open file, #'//trim(string(ios)))
!detect FE topology and count number of elements and nodes
nnod = 0
do
  read (unit=iu, fmt='(a8)', iostat=ios) str
  if (ios == iostat_end) then
    exit !end of file
  elseif (ios /= 0) then
    call error('(module_patran/load_patran) unable to read from file, #'//trim(string(ios)))
  end if
  do i = 1, size(t,1)
    if (str == t(i)%name) then !elements
      if (len_trim(elto) == 0 .or. fedim < t(i)%fedim) then !never saved or bigger dimension than saved
        elto = t(i)%name
        fedim = t(i)%fedim
        itopo = t(i)%itopo
        !dim is stablished later
        lnn   = t(i)%lnn !npe
        lnv   = t(i)%lnv !vpe
        lne   = t(i)%lne !ape
        lnf   = t(i)%lnf !cpe
        nel   = 1
        exit
      elseif (elto == t(i)%name) then !already saved
        nel = nel + 1
        exit
      else
        call error('(module_patran/load_patran) program cannot deal with 2 FE type of the same dimension: '//&
        &trim(elto)//' and '//trim(t(i)%name))
      end if
    elseif (str == 'GRID') then !nodes
      nnod = nnod + 1
    end if
  end do
enddo
if (fedim == 0) call error('(module_patran/load_patran) FE topolgy not found')
if (nnod  == 0) call error('(module_patran/load_patran) nodes not found')
rewind(iu)
!other variables
nmat = 1
elto = ' '

 !! mm, MATRIZ DE CONECTIVIDADES Y nsd, NUMERO DE SUBDOMINIO
    allocate(nsd(nel))
    select case(itopo)
      case(1)
        vpe=3; npe=3; ape=3; cpe=1; topo='CTRIA3'; nver=nnod
        allocate(mm(npe,nel))

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,201)elto,h,nsd(h),(mm(k,h),k=1,npe)
            read(iu,100)elto
            backspace(iu)
          end do
        end do


      case(2)
        vpe=3; npe=6; ape=3; cpe=1; topo='CTRIA6'
        allocate(mmver(nnod),nn(npe,nel),mm(3,nel))

        nver=0
        do j=1,nnod
          mmver(j)=0
        end do

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,203)elto,h,nsd(h),(nn(k,h),k=1,npe)
            read(iu,100)elto
            backspace(iu)

            do k=1,3
              if (mmver(nn(k,h)).eq.0) then
                nver=nver+1
                mmver(nn(k,h))=nver
                mm(k,h)=nver
              else
                mm(k,h)=mmver(nn(k,h))
              end if
            end do

          end do
        end do


      case(3)
        vpe=3; npe=3; ape=3; cpe=1; topo='CTRIA6'
        allocate(mm(3,nel),mmver(nnod),mmnod(nnod),nn(3,nel),m(6,nel))

        nver=0;nnnod=0
        do j=1,nnod
          mmver(j)=0
          mmnod(j)=0
        end do

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,203)elto,h,nsd(h),(m(k,h),k=1,6)
            read(iu,100)elto
            backspace(iu)

            do k=1,3
              if (mmver(m(k,h)).eq.0) then
                nver=nver+1
                mmver(m(k,h))=nver
                mm(k,h)=nver
              else
                mm(k,h)=mmver(m(k,h))
              end if
              kk=k+3
              if (mmnod(m(kk,h)).eq.0) then
                nnnod=nnnod+1
                mmnod(m(kk,h))=nnnod
                nn(k,h)=nnnod
              else
                nn(k,h)=mmnod(m(kk,h))
              end if
            end do

          end do

        end do


      case(4)
        vpe=4; npe=4; ape=4; cpe=1; topo='CQUAD4'; nver=nnod
        allocate(mm(npe,nel))

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,202)elto,h,nsd(h),(mm(k,h),k=1,npe)
            read(iu,100)elto
            backspace(iu)
          end do
        end do


      case(5)
        vpe=4; npe=8; ape=4; cpe=1; topo='CQUAD8'
        allocate(mmver(nnod),nn(npe,nel),mm(4,nel))

        nver=0
        do j=1,nnod
          mmver(j)=0
        end do

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,204)elto,h,nsd(h),(nn(k,h),k=1,npe)
            read(iu,212)(nn(k,h),k=7,8)
            read(iu,100)elto
            backspace(iu)

            do k=1,4
              if (mmver(nn(k,h)).eq.0) then
                nver=nver+1
                mmver(nn(k,h))=nver
                mm(k,h)=nver
              else
                mm(k,h)=mmver(nn(k,h))
              end if
            end do

          end do
        end do


      case(6)
        vpe=4; npe=4; ape=4; cpe=1; topo='CQUAD8'
        allocate(mm(4,nel),mmver(nnod),mmnod(nnod),nn(4,nel),m(8,nel))

        nver=0;nnnod=0
        do j=1,nnod
          mmver(j)=0
          mmnod(j)=0
        end do

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,204)elto,h,nsd(h),(m(k,h),k=1,8)
            read(iu,100)elto
            backspace(iu)

            do k=1,4
              if (mmver(m(k,h)).eq.0) then
                nver=nver+1
                mmver(m(k,h))=nver
                mm(k,h)=nver
              else
                mm(k,h)=mmver(m(k,h))
              end if
              kk=k+4
              if (mmnod(m(kk,h)).eq.0) then
                nnnod=nnnod+1
                mmnod(m(kk,h))=nnnod
                nn(k,h)=nnnod
              else
                nn(k,h)=mmnod(m(kk,h))
              end if
            end do

          end do
        end do


      case(7)
        vpe=4; npe=4; ape=6; cpe=4; topo='CTETRA'; nver=nnod
        allocate(mm(npe,nel))

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)
          do while(elto.eq.topo)
            read(iu,202)elto,h,nsd(h),(mm(k,h),k=1,npe)
            read(iu,100)elto
            backspace(iu)
          end do
        end do


      case(8)
        vpe=6; npe=6; ape=9; cpe=5; topo='CPENTA'; nver=nnod
        allocate(mm(npe,nel))

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,203)elto,h,nsd(h),(mm(k,h),k=1,npe)
            read(iu,100)elto
            backspace(iu)
          end do
        end do


      case(9)
        vpe=8; npe=8; ape=12; cpe=6; topo='CHEXA'; nver=nnod
        allocate(mm(npe,nel))

        do i=1,nmat
          do while(elto.ne.topo)
            read(iu,100)elto
          end do
          backspace(iu)

          do while(elto.eq.topo)
            read(iu,203)elto,h,nsd(h),(mm(k,h),k=1,6)
            read(iu,212)(mm(k,h),k=7,8)
            read(iu,100)elto
            backspace(iu)
          end do
        end do

    end select

!Modificacion_Fran
rewind(iu)
!Fin_Modificacion_Fran

 !! z, COORDENADAS

    grid=' '
    read(iu,99)grid

    do while(grid.ne.'GRID')
      read(iu,99)grid
    end do
    backspace(iu)

    if ((itopo.eq.1).or.(itopo.eq.4)) then
      allocate(z(2,nnod))
      do while(grid.eq.'GRID')
        read(iu,100)gridd
        backspace(iu)
        if (gridd.eq.'GRID') then
          read(iu,221)gridd,i,(z(k,i),k=1,2)
          read(iu,99)grid
        else
          read(iu,224)gridd,i,(z(k,i),k=1,2)
          read(iu,*)
          read(iu,99)grid
        end if
        backspace(iu)
      end do
    end if

    if ((itopo.eq.2).or.(itopo.eq.5)) then
      allocate(z(2,nnod))
      do while(grid.eq.'GRID')
        read(iu,100)gridd
        backspace(iu)
        if (gridd.eq.'GRID') then
          read(iu,221)gridd,i,(z(k,mmver(i)),k=1,2)
          read(iu,99)grid
        else
          read(iu,224)gridd,i,(z(k,mmver(i)),k=1,2)
          read(iu,*)
          read(iu,99)grid
        end if
        backspace(iu)
      end do
    end if

    if ((itopo.eq.3).or.(itopo.eq.6)) then
      allocate(zn(2,nnod),z(2,nver),zz(2,nnnod))
      do while(grid.eq.'GRID')
        read(iu,100)gridd
        backspace(iu)
        if (gridd.eq.'GRID') then
          read(iu,221)gridd,i,(zn(k,i),k=1,2)
          read(iu,99)grid
        else
          read(iu,224)gridd,i,(zn(k,i),k=1,2)
          read(iu,*)
          read(iu,99)grid
        end if
        backspace(iu)
      end do
      do k=1,nnod
        if (mmver(k).ne.0) then
          z(1,mmver(k))=zn(1,k)
          z(2,mmver(k))=zn(2,k)
        end if

        if (mmnod(k).ne.0) then
          zz(1,mmnod(k))=zn(1,k)
          zz(2,mmnod(k))=zn(2,k)
        end if
      end do
    end if

    if (itopo.gt.6) then
      allocate(z(3,nnod))
      do while(grid.eq.'GRID')
        read(iu,100)gridd
        backspace(iu)
        if (gridd.eq.'GRID') then
          read(iu,222)gridd,i,(z(k,i),k=1,3)
          read(iu,99)grid
        else
          read(iu,224)gridd,i,(z(k,i),k=1,2)
          read(iu,223)z(3,i)
          read(iu,99)grid
        end if
        backspace(iu)
      end do
    end if

    !Modificacion_Fran
    !si det < 0, se intercambian los dos ultimos nodos del elemento
    DIM = 2; if (itopo > 6) DIM = 3
    if (itopo == 1 .or. itopo == 7) then
      do k = 1, nel
        if (det(reshape( (/ ((z(j,mm(i,k))-z(j,mm(1,k)), j=1,DIM), i=2,DIM+1) /), &
        (/DIM, DIM/) )) < 0 ) then
          tmp = mm(vpe-1, k); mm(vpe-1, k) = mm(vpe, k); mm(vpe, k) = tmp
          if (det(reshape( (/ ((z(j,mm(i,k))-z(j,mm(1,k)), j=1,DIM), i=2,DIM+1) /), &
          (/DIM, DIM/) )) < 0 ) stop 'El programa obtiene determinantes negativos'
        end if
      end do
    end if
    !Fin_Modificacion_Fran

 !! nrv, NUMERO DE REFERENCIA DE VERTICES

    allocate(nnrv(nnod),nrv(vpe,nel))
    do j=1,nnod
      nnrv(j)=0
    end do

    hd=0

    !Interpretaci�n de los SPC
    spc=' '
    read(iu,100)spc

    !Modificacion_Fran
    !do while(spc.ne.'SPC1')
    do while(spc.ne.'SPC1' .and. trim(spc)/='SPC')
    !Fin_Modificacion_Fran
      read(iu,100)spc
      if (spc.eq.'ENDDATA') then
        rewind(iu)
        exit
      end if
    end do
    backspace(iu)

    !Modificacion_Fran
    !do while((spc.eq.'SPC1').or.(spc.eq.'   '))
    do while((spc.eq.'SPC1').or.(spc.eq.'   ').or. trim(spc) == 'SPC')
    !Fin_Modificacion_Fran
      if (spc.eq.'SPC1') then
        read(iu,204,iostat=status)spc,h,g,(nodo(k),k=1,8)
        if (h.gt.hd) then
          hd=h
        end if
      !Modificacion_Fran
      elseif (trim(spc) == 'SPC') then
        call set_SPC(iu)
      !Fin_Modificacion_Fran
      else
        read(iu,215,iostat=status)(nodo(k),k=1,8)
      end if

      if (status.gt.0) then   ! Hay un THROUGH
        backspace(iu)
        read(iu,225)spc,h,g,nodo(1),through,nodo(3)

        do k=nodo(1),nodo(3)
          nnrv(k)=h
        end do
      !Modificacion_Fran
      elseif (trim(spc) == 'SPC') then
        !la asignacion se realiza una vez leidas todas las condiciones,
        !cuando se sepa agrupar los nodos que corresponden a la misma condicion
        continue
      !Fin_Modificacion_Fran
      else
        k=2
        do while (k.lt.9)
          if (nodo(k).eq.0) then
            exit
          end if
          k=k+1
        end do
        if (k.lt.9) then
          select case(k)
            case(2)
              !Modificacion_Fran
              !nnrv(nodo(1))=h
              nnrv(nodo(1)) = max(nnrv(nodo(1)), h)
              !Fin_Modificacion_Fran
            case(3)
              !Modificacion_Fran
              !nnrv(nodo(1))=h
              !nnrv(nodo(2))=h
              nnrv(nodo(1)) = max(nnrv(nodo(1)), h)
              nnrv(nodo(2)) = max(nnrv(nodo(2)), h)
              !Fin_Modificacion_Fran
            case(4)
              do k=1,3
                !Modificacion_Fran
                !nnrv(nodo(k))=h
                nnrv(nodo(k)) = max(nnrv(nodo(k)), h)
                !Fin_Modificacion_Fran
              end do
            case(5)
              do k=1,4
                !Modificacion_Fran
                !nnrv(nodo(k))=h
                nnrv(nodo(k)) = max(nnrv(nodo(k)), h)
                !Fin_Modificacion_Fran
              end do
            case(6)
              do k=1,5
                !Modificacion_Fran
                !nnrv(nodo(k))=h
                nnrv(nodo(k)) = max(nnrv(nodo(k)), h)
                !Fin_Modificacion_Fran
              end do
            case(7)
              do k=1,6
                !Modificacion_Fran
                !nnrv(nodo(k))=h
                nnrv(nodo(k)) = max(nnrv(nodo(k)), h)
                !Fin_Modificacion_Fran
              end do
            case(8)
              do k=1,7
                !Modificacion_Fran
                !nnrv(nodo(k))=h
                nnrv(nodo(k)) = max(nnrv(nodo(k)), h)
                !Fin_Modificacion_Fran
              end do
          end select
        else
          do k=1,8
            !Modificacion_Fran
            !nnrv(nodo(k))=h
            nnrv(nodo(k)) = max(nnrv(nodo(k)), h)
            !Fin_Modificacion_Fran
          end do
        end if
      end if

      read(iu,100)spc
      !Modificacion_Fran
      !if ((spc.ne.'SPC1').and.(spc.ne.'  ')) then
      if ((spc.ne.'SPC1').and.(spc.ne.'  ') .and. trim(spc) /= 'SPC') then
      !Fin_Modificacion_Fran
        read(iu,100)spc
      end if
      backspace(iu)
    end do
    !Modificacion_Fran
    call assign_SPC(nnrv, hd)
    hd = maxval(nnrv)
    !Fin_Modificacion_Fran

    !Interpretaci�n de los FORCE
    force=' '
    backspace(iu)
    read(iu,100)force

    do while(force.ne.'FORCE')
      read(iu,100)force

      if (force.eq.'ENDDATA') then
        exit
      end if

    end do
    backspace(iu)

    do while(force.eq.'FORCE')

      !Modificacion_Fran
      !read(iu,212)h,g
      !if ((nnrv(g).eq.0).or.(nnrv(g).gt.hd)) then
      !nnrv(g)=hd+h
      !end if
      call set_FORCE(iu)
      !Fin_Modificacion_Fran

      read(iu,100)force
      if (force.ne.'FORCE') then
        read(iu,100)force
      end if
      backspace(iu)
      if (force.eq.'MOMENT') then
        do while (force.eq.'MOMENT')
          read(iu,100)force
        end do
        read(iu,100)force
        backspace(iu)
      end if
    end do
    !Modificacion_Fran
    call assign_FORCE(nnrv, hd)
    !call save(trim(fileout)//'.rff', 12, nver, real(nnrv))
    !Fin_Modificacion_Fran


    !nrv
    select case(itopo)
      case(2,5)
        do j=1,nel
          do k=1,vpe
            nrv(k,j)=nnrv(nn(k,j))
          end do
        end do

      case(3,6)
        do j=1,nel
          do k=1,vpe
            nrv(k,j)=nnrv(m(k,j))
          end do
        end do

      case default
        do j=1,nel
          do k=1,vpe
            nrv(k,j)=nnrv(mm(k,j))
          end do
        end do
    end select


 !! nra, NUMERO DE REFERENCIA DE ARISTAS

    allocate(nra(ape,nel))

   ! DOS DIMENSIONES

    if (itopo.lt.7) then
      h=0   !contador de aristas vistas
      do j=1,nel
        do k=1,ape-1
          nra(k,j)=dos(nrv(k,j),nrv(k+1,j),hd)
        end do
        nra(ape,j)=dos(nrv(ape,j),nrv(1,j),hd)
      end do

!**********************************************************************************************************************************
      !Modificacion_Fran: lectura de nra mediante CBEAM, para CTRIA3
      !Cada arista CBEAM tiene un número de referencia asociado, que debemos guardar en nra
      !Se usa el mecanismo creado para unv2mfm para pasar la información a nrc
      if (itopo==1) then !los EF bidimensionales son triangulos Lagrange P1
        rewind(iu) !se buscan los elementos desde el principio del fichero
        call alloc(groups, 1) !alojaremos los elementos CBEAM en groups(1)
        do while (.true.) !recoleccion de elementos CBEAM
          read(iu,fmt=*, iostat=io) typef
          if (io == iostat_end) exit
          if (trim(adjustl(typef)) == 'CBEAM') then
            backspace(iu, iostat=io); if (io /= 0) exit
            read(iu,fmt='(a8,4i8)', iostat=io) typef, ind, ref, (vnum(i), i=1, 2); if (io /= 0) exit
            call insert(groups(1), [order(vnum(1:2)), ref])
           endif
        enddo
        if (groups(1)%ncells > 0) then !construccion de nrc
          do k = 1, nel !bucle en tetraedros
            do j = 1, 3 !bucle en aristas
              pos = search(groups(1), order(mm(edge(:,j),k)))
              if (pos > 0) nra(j,k) = groups(1)%cell(pos,3)
            enddo
          enddo
        endif
      endif
      !Modificacion_Fran: lectura de nra mediante CBEAM, para CTRIA6
      !Cada arista CBEAM tiene un número de referencia asociado, que debemos guardar en nra
      !Se usa el mecanismo creado para unv2mfm para pasar la información a nrc
      if (itopo==2) then !los EF bidimensionales son triangulos Lagrange P2
        rewind(iu) !se buscan los elementos desde el principio del fichero
        call alloc(groups, 1) !alojaremos los elementos CBEAM en groups(1)
        do while (.true.) !recoleccion de elementos CBEAM
          read(iu,fmt=*, iostat=io) typef
          if (io == iostat_end) exit
          if (trim(adjustl(typef)) == 'CBEAM') then
            backspace(iu, iostat=io); if (io /= 0) exit
            read(iu,fmt='(a8,4i8)', iostat=io) typef, ind, ref, (vnum(i), i=1, 2); if (io /= 0) exit
            call insert(groups(1), [order(vnum(1:2)), ref])
           endif
        enddo
        if (groups(1)%ncells > 0) then !construccion de nrc
          do k = 1, nel !bucle en tetraedros
            do j = 1, 6 !bucle en aristas
              pos = search(groups(1), order(nn(edge2(:,j),k)))
              if (pos > 0) nra((j-1)/2+1,k) = groups(1)%cell(pos,3)
            enddo
          enddo
        endif
      endif
      !Fin_Modificacion_Fran
!**********************************************************************************************************************************
    end if


   ! TETRAEDROS

    if (itopo.eq.7) then
      h=0
      do j=1,nel
        do k=1,2
          nra(k,j)=dos(nrv(k,j),nrv(k+1,j),hd)
        end do
        nra(3,j)=dos(nrv(3,j),nrv(1,j),hd)

        do k=1,3
          nra(3+k,j)=dos(nrv(k,j),nrv(4,j),hd)
        end do
      end do
    end if

    if (itopo.eq.8) then
      do j=1,nel
        do k=1,2
          nra(k,j)=dos(nrv(k,j),nrv(k+1,j),hd)
        end do
        nra(3,j)=dos(nrv(3,j),nrv(1,j),hd)

        do k=1,3
          nra(3+k,j)=dos(nrv(k,j),nrv(3+k,j),hd)
        end do

        do k=4,5
          nra(3+k,j)=dos(nrv(k,j),nrv(k+1,j),hd)
        end do
        nra(9,j)=dos(nrv(6,j),nrv(4,j),hd)
      end do
    end if

    if (itopo.eq.9) then
      do j=1,nel
        do k=1,3
          nra(k,j)=dos(nrv(k,j),nrv(k+1,j),hd)
        end do
        nra(4,j)=dos(nrv(4,j),nrv(1,j),hd)

        do k=1,4
          nra(4+k,j)=dos(nrv(k,j),nrv(4+k,j),hd)
        end do

        do k=5,7
          nra(4+k,j)=dos(nrv(k,j),nrv(k+1,j),hd)
        end do
        nra(12,j)=dos(nrv(8,j),nrv(5,j),hd)
      end do
    end if



 !! nrc, N�MERO DE REFERENCIA DE CARAS

   allocate(nrc(cpe,nel))
   nrc = 0

   !tetraedro
    if (itopo.eq.7) then
      h=0
      !Modificacion_Fran
      !do j=1,nel
      !  nrc(1,j)=tres(nrv(1,j),nrv(2,j),nrv(3,j),hd)
      !  nrc(2,j)=tres(nrv(1,j),nrv(3,j),nrv(4,j),hd)
      !  nrc(3,j)=tres(nrv(1,j),nrv(2,j),nrv(4,j),hd)
      !  nrc(4,j)=tres(nrv(2,j),nrv(3,j),nrv(4,j),hd)
      !end do

      !if (caras_interiores /= 's') call erase_interior_faces(nel, mm, nrc)

    !Modificacion_Fran
    !Modificacion para la lectura de nrc en caso de tetraedros Tet4
    !Se supone que tenemos tetraedros y que sus caras son CTRIA3. Cada
    !cara lleva un número de referncia asociado, que debemos guardar en nrc
    !se usa el mecanismo creado para unv2mfm para pasar la información a nrc

      if (itopo==7) then !los EF tridimensionales son tetraedros Lagrange P1
        rewind(unit=iu) !se buscan los elementos desde el principio del fichero
        call alloc(groups, 1) !alojaremos los elementos CTRIA3 en groups(1)
        do while (.true.) !recoleccion de elementos CTRIA3
          read(unit=iu,fmt=*, iostat=io) typef
          if (io == iostat_end) exit
          if (trim(adjustl(typef)) == 'CTRIA3') then
            backspace(iu, iostat=io)
            if (io /= 0) exit
            read(unit=iu,fmt=*, iostat=io) typef, ind, ref, (vnum(i), i=1, 3)
            if (io /= 0) exit
            call insert(groups(1), [order(vnum), ref])
           endif
        enddo
        if (groups(1)%ncells > 0) then !construccion de nrc
          do k = 1, nel !bucle en tetraedros
            do j = 1, 4 !bucle en caras
              pos = search(groups(1), order(mm(face(:,j),k)))
              if (pos > 0) nrc(j,k) = groups(1)%cell(pos,4)
            enddo
          enddo
        endif
      endif
      !Fin_Modificacion_Fran
    end if

   !prisma
    if (itopo.eq.8) then
      h=0
      do j=1,nel
        nrc(1,j)=tres(nrv(1,j),nrv(2,j),nrv(3,j),hd)
        nrc(2,j)=cuatro(nrv(1,j),nrv(3,j),nrv(4,j),nrv(6,j),hd)
        nrc(3,j)=cuatro(nrv(1,j),nrv(2,j),nrv(4,j),nrv(5,j),hd)
        nrc(4,j)=tres(nrv(4,j),nrv(5,j),nrv(6,j),hd)
        nrc(5,j)=cuatro(nrv(2,j),nrv(3,j),nrv(5,j),nrv(6,j),hd)
      end do
    end if

   !hexaedro
    if (itopo.eq.9) then
      h=0
      do j=1,nel
        nrc(1,j)=cuatro(nrv(1,j),nrv(2,j),nrv(3,j),nrv(4,j),hd)
        nrc(2,j)=cuatro(nrv(1,j),nrv(4,j),nrv(5,j),nrv(8,j),hd)
        nrc(3,j)=cuatro(nrv(1,j),nrv(2,j),nrv(5,j),nrv(6,j),hd)
        nrc(4,j)=cuatro(nrv(5,j),nrv(6,j),nrv(7,j),nrv(8,j),hd)
        nrc(5,j)=cuatro(nrv(2,j),nrv(3,j),nrv(6,j),nrv(7,j),hd)
        nrc(6,j)=cuatro(nrv(4,j),nrv(3,j),nrv(7,j),nrv(8,j),hd)
      end do
    end if

  close(iu)

print'(a,i9)','File loaded!'
print'(a,i9)','Global number of elements: ', nel
print'(a,i9)','Global number of nodes:    ', nnod
print'(a,i9)','Global number of vertices: ', nver
print'(a,i9)','Space dimension :          ', dim
print'(a,i9)','Local number of nodes :    ', lnn
print'(a,i9)','Local number of vertices : ', lnv
print'(a,i9)','Local number of edges :    ', lne
print'(a,i9)','Local number of faces :    ', lnf


!! ARCHIVO DE SALIDA
!  if (form.eq.1) then
!    if (itopo.lt.7) then
!      if ((itopo.eq.1).or.(itopo.eq.4)) then
!        open(4,file=fileout,form='formatted')
!          rewind(4)
!          write(4,*)nel,nnod,nver
!          write(4,*)((mm(k,j),k=1,npe),j=1,nel), &
!            &  ((nra(k,j),k=1,ape),j=1,nel), ((nrv(k,j),k=1,vpe),j=1,nel), &
!            &  ((z(k,j),k=1,2),j=1,nver)
!          write(4,*)(nsd(j),j=1,nel)
!        close(4)
!      end if

!      if ((itopo.eq.2).or.(itopo.eq.5)) then
!        open(4,file=fileout,form='formatted')
!          rewind(4)
!          write(4,'(8(i11))')nel,nnod,nver,2,npe,vpe,ape,1
!          write(4,*)((nn(k,j),k=1,npe),j=1,nel),((mm(k,j),k=1,vpe),j=1,nel), &
!            &  ((nra(k,j),k=1,ape),j=1,nel), ((nrv(k,j),k=1,vpe),j=1,nel), &
!            &  ((z(k,j),k=1,2),j=1,nver)
!          write(4,*)(nsd(j),j=1,nel)
!        close(4)
!      end if

!      if ((itopo.eq.3).or.(itopo.eq.6)) then
!        open(4,file=fileout,form='formatted')
!          rewind(4)
!          write(4,*)nel,nnnod,nver
!          write(4,*)((nn(k,j),k=1,npe),j=1,nel),((mm(k,j),k=1,npe),j=1,nel), &
!            &  ((nra(k,j),k=1,ape),j=1,nel), ((nrv(k,j),k=1,vpe),j=1,nel), &
!            &  ((z(k,j),k=1,2),j=1,nver)
!          write(4,*)(nsd(j),j=1,nel)
!          write(4,*)((zz(k,j),k=1,2),j=1,nnnod)
!        close(4)
!      end if

!    else !casos 7, 8 y 9
!      open(4,file=fileout,form='formatted')
!        rewind(4)
!        write(4,'(8(i11))')nel,nnod,nver,3,npe,npe,ape,cpe
!        write(4,*)((mm(k,j),k=1,npe),j=1,nel),((nrc(k,j),k=1,cpe),j=1,nel), &
!          &  ((nra(k,j),k=1,ape),j=1,nel),((nrv(k,j),k=1,vpe),j=1,nel), &
!          &  ((z(k,j),k=1,3),j=1,nver)
!        write(4,*)(nsd(j),j=1,nel)
!      close(4)
!    end if

!  else
!    if (itopo.lt.7) then
!      if ((itopo.eq.1).or.(itopo.eq.4)) then
!        open(4,file=fileout,form='unformatted')
!          rewind(4)
!          write(4)nel,nnod,nver
!          write(4)((mm(k,j),k=1,npe),j=1,nel), &
!            &  ((nra(k,j),k=1,ape),j=1,nel), ((nrv(k,j),k=1,vpe),j=1,nel), &
!            &  ((z(k,j),k=1,2),j=1,nver)
!          write(4)(nsd(j),j=1,nel)
!        close(4)
!      end if

!      if ((itopo.eq.2).or.(itopo.eq.5)) then
!        open(4,file=fileout,form='unformatted')
!          rewind(4)
!          write(4)nel,nnod,nver
!          write(4)((nn(k,j),k=1,npe),j=1,nel),((mm(k,j),k=1,vpe),j=1,nel), &
!            &  ((nra(k,j),k=1,ape),j=1,nel), ((nrv(k,j),k=1,vpe),j=1,nel), &
!            &  ((z(k,j),k=1,2),j=1,nver)
!          write(4)(nsd(j),j=1,nel)
!        close(4)
!      end if

!      if ((itopo.eq.3).or.(itopo.eq.6)) then
!        open(4,file=fileout,form='unformatted')
!          rewind(4)
!          write(4)nel,nnnod,nver
!          write(4)((nn(k,j),k=1,npe),j=1,nel),((mm(k,j),k=1,npe),j=1,nel), &
!            &  ((nra(k,j),k=1,ape),j=1,nel), ((nrv(k,j),k=1,vpe),j=1,nel), &
!            &  ((z(k,j),k=1,2),j=1,nver)
!          write(4)(nsd(j),j=1,nel)
!          write(4)((zz(k,j),k=1,2),j=1,nnnod)
!        close(4)
!      end if

!    else ! casos 7, 8 y 9
!      open(4,file=fileout,form='unformatted')
!        rewind(4)
!        write(4)nel,nnod,nver
!        write(4)((mm(k,j),k=1,npe),j=1,nel),((nrc(k,j),k=1,cpe),j=1,nel),   &
!          &  ((nra(k,j),k=1,ape),j=1,nel),((nrv(k,j),k=1,vpe),j=1,nel), &
!          &  ((z(k,j),k=1,3),j=1,nver)
!        write(4)(nsd(j),j=1,nel)
!      close(4)
!    end if
!  end if




 99  format(a4)
100  format(a8)

201 format(a8,5i8)
202 format(a8,6i8)
203 format(a8,8i8)
204 format(a8,10i8)

!211 format(8x,i8)
212 format(8x,2i8)
!213 format(8x,4i8)
!214 format(8x,6i8)
215 format(8x,8i8)

221 format(a8,i8,8x,2d8.5)
222 format(a8,i8,8x,3d8.5)
223 format(8x,d16.5)
224 format(a8,i16,16x,d16.5,d16.5)
225 format(a8,3i8,a8,i8)

end subroutine

!Modificacion_Fran

!-----------------------------------------------------------------------
! order: order an array (DIM <=3)
!-----------------------------------------------------------------------
function order(u) result(v)
integer, dimension(:), intent(in) :: u
integer, dimension(size(u,1)) :: v

v = u
if (size(u,1)==2 .and. v(1) > v(2)) then
  call swap(v(1), v(2))
elseif (size(u,1)==3) then
  if (v(1) > v(2)) call swap(v(1), v(2))
  if (v(2) > v(3)) call swap(v(2), v(3))
  if (v(1) > v(2)) call swap(v(1), v(2))
end if

end function

!-----------------------------------------------------------------------
! swap: swap two variables
!-----------------------------------------------------------------------
subroutine swap(u,v)
integer, intent(inout) :: u, v
integer :: tmp

tmp = u; u = v; v = tmp

end subroutine
!Fin_Modificacion_Fran

end module
