module module_ansys
!-----------------------------------------------------------------------
! Utility to convert ANSYS mesh files into MFM Lagrange P1 
! tetrahedra mesh
!
! Licensing: This code is distributed under the GNU GPL license.
! Authors: Iban Constenla
!          Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 24/05/2013
!
! PUBLIC PROCEDURES:
!   load_ansys: loads an ANSYS mesh file
!
! PRIVATE PROCEDURES:
!   load_ansys_tria: loads an ANSYS triangular mesh file
!   load_ansys_tria: loads an ANSYS tetrahedral mesh file
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64, iostat_end
use module_os_dependant, only: maxpath
use module_report, only: error
use module_convers, only: string, word_count, word, int
use module_alloc_int_r1, only: add, reduce
use module_set, only: sunique
implicit none

!Private procedures
private :: load_ansys_tria, load_ansys_tetra

contains

!-----------------------------------------------------------------------
! load_ansys: loads an ANSYS tetrahedra mesh file
!-----------------------------------------------------------------------
subroutine load_ansys(file_in, iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
character(*),              intent(in)  :: file_in
integer,                   intent(in)  :: iu
integer,                   intent(out) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
integer,      allocatable, intent(out) :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64), allocatable, intent(out) :: z(:,:)
integer :: ios, nets
integer, allocatable :: et(:), types(:)
character(maxpath) :: str, eltype

!open file
open (unit=iu, file=file_in, form='formatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('(module_ansys/load_ansys) unable to open file, #'//trim(string(ios)))
!store element types from sections n. 12
nets = 0
do
  read (unit=iu, fmt='(a)', iostat=ios) str
  if (ios == iostat_end) then
    exit !end of file
  elseif (ios /= 0) then
    call error('(module_ansys/load_ansys) unable to read from file, #'//trim(string(ios)))
  end if
  if (word_count(str) >= 2) then
    if (trim(word(str,1)) == '(12' .and. trim(word(str,2)) /= '(0') then
      eltype = word(str,6)
      nets = nets+1
      call add(et, int(eltype(1:len_trim(eltype)-2)), nets, fit=.false.)
    end if
  end if
end do
call reduce(et, nets)
!check whether there is a single element type
call sunique(et, types)
if (size(types,1) /= 1) call error('(module_ansys/load_ansys) file stores more than an element type')

select case (types(1))
case(1) !triangular
  call load_ansys_tria(iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nra, nrv, z, nsd)
case(2) !tetrahedra
  call load_ansys_tetra(iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
case default
  call error('(module_ansys/load_ansys) ANSYS element type '//trim(string(types(1)))//' not implement')
end select

end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! load_ansys_tria: loads an ANSYS triangular mesh file
!-----------------------------------------------------------------------
subroutine load_ansys_tria(iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nra, nrv, z, nsd)
integer,                   intent(in)  :: iu
integer,                   intent(out) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
integer,      allocatable, intent(out) :: nn(:,:), mm(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64), allocatable, intent(out) :: z(:,:)

integer, dimension(:),   allocatable :: contador,bref,nrvg,lfn
integer, dimension(:,:), allocatable :: cver,celad
integer :: ndcar, nref, nsref, nnoel, npuel, narel, ncael, fi, li, i, j, k, k_car, idx, &
st, bctype, ncarint, refname, temp, ios
real(real64), allocatable :: det(:)
real(real64) :: x(4), y(4)
character(60) :: idvar, var, hex_ver, hex_car, hex_elem, f_i, l_i, hextype, len
character(20) :: aux(5)
character(maxpath) :: str, wd

rewind(iu)
! Comments 
  do i=1,7
    read(iu,'(a)')
  end do
! Global dimensions
read (unit=iu, fmt='(a)', iostat=ios) str
if (ios /= 0) call error('(module_ansys/load_ansys_tria) unable to read from file, #'//trim(string(ios)))
if (word_count(str) < 2) call error('(module_ansys/load_ansys_tria) unable to read idx, dim')
wd = word(str,1); idx = int(wd(2:len_trim(wd)))
wd = word(str,2); dim = int(wd(1:len_trim(wd)-1))

! Hexadecimal reading
  read(iu,*) var,var,var,hex_ver,var
  read(iu,*) var,var,var,hex_car,var
  read(iu,*) var,var,var,hex_elem,var
! Global vertices
  read(hex_ver,'(z20)') nver
! Global faces
  read(hex_car,'(z20)') ndcar
! Global elements
  read(hex_elem,'(z20)') nel  
  write(len,*) 1+int(log(float(nel))/log(16.))
  len='z'//trim(adjustl(len))  
! Allocatables
  allocate(z(dim,nver),cver(dim,ndcar),celad(2,ndcar))
  allocate(bref(ndcar),nsd(nel))
! Initialization
  bref(1:ndcar)=0
  cver(1:dim,1:ndcar)=0
  celad(1:2,1:ndcar)=0
  nref=1  
  nsref=1
  refname=1
! 2D meshes
!print*,'========================================'
!print*,'References :'
!print*,'========================================'
! Loop to EOF (end of file)
bc: do while (.true.)
      read(iu,*,iostat=st) idvar
      if (is_iostat_end(st)) exit
          if (trim(idvar) /= '(45') then
             backspace(iu)
             read(iu,*) idvar,var,f_i,l_i
             read(f_i,'(z20)') fi
             read(l_i,'(z20)') li
             if (trim(idvar) == '(10') then
                ! Coordinates of vertices
                do j=fi,li
                  read(iu,*) z(1:dim,j)
                end do
                read(iu,'(a)')          
             elseif (trim(idvar) == '(13') then
                backspace(iu)
                read(iu,*) idvar,var,f_i,l_i,hextype
                read(hextype,'(z20)') bctype
                ! Inside faces
                if (bctype == 2) then
                    do j=fi,li
                       read(iu,*) aux(1:4)
                       read(aux(1),'(z20)') cver(1,j)
                       read(aux(2),'(z20)') cver(2,j)
                       read(aux(3),'(z20)') celad(1,j)
                       read(aux(4),'(z20)') celad(2,j)
                       ncarint=li
                    end do
                ! Frontier references                
                else
                    do j=fi,li
                       read(iu,*) aux(1:4)
                       read(aux(1),'(z20)') cver(1,j)
                       read(aux(2),'(z20)') cver(2,j)
                       read(aux(3),'(z20)') celad(1,j)
                       read(aux(4),'(z20)') celad(2,j)
                       bref(j)=nref
                    end do
                    nref=nref+1
                end if
                read(iu,'(a)')
                read(iu,'(a)')
             elseif (trim(idvar) == '(12') then
                ! Subdomains
                do j=fi,li
                  nsd(j)=nsref
                end do
                nsref=nsref+1
             end if
          else
             backspace(iu)
             read(iu,*) idvar,var,f_i,l_i
             if (trim(f_i) /= 'interior' .and. trim(f_i) /= 'fluid') then
              ! Print references
!                print*,l_i(1:len_trim(l_i)-4),' Correspond to reference number :',refname
               print'(a)','Reference number '//trim(string(refname))//', associated with group named "'//l_i(1:len_trim(l_i)-4)//'"'
                refname=refname+1
             else  
                cycle bc
             end if 
          end if
        end do bc
      close(iu)
!      print*,'========================================'  
!      print*,'Calculating ... '
!      print*,'========================================'  
!**************************************************************************
! mm: integer matrix (3,nel) where mm(i,k) is the global index of the i-th
! node of the k-th element of the mesh
!**************************************************************************
  allocate(mm(3,nel),contador(nel))
  mm(1:3,1:nel)=0
  contador(1:nel)=0
  do j=1,ncarint
    do i=1,2
      k_car=celad(i,j)
      if (contador(k_car)==0) then
        mm(1:2,k_car)=cver(1:2,j)        
        contador(k_car)=1
      else if (contador(k_car)==1) then
        do k=1,2
          if (cver(k,j) /= mm(1,k_car) .and. &
              cver(k,j) /= mm(2,k_car)) then
           mm(3,k_car) = cver(k,j)
           contador(k_car)=2
           exit
          end if
        end do
      end if
    end do      
  end do
  do j=ncarint+1,ndcar
    k_car=celad(1,j)
    if (contador(k_car)==0) then
      mm(1:2,k_car)=cver(1:2,j)        
      contador(k_car)=1
    else if (contador(k_car)==1) then
      do k=1,2
        if (cver(k,j) /= mm(1,k_car) .and. &
            cver(k,j) /= mm(2,k_car)) then
         mm(3,k_car) = cver(k,j)
         contador(k_car)=2
         exit
        end if
      end do
    end if    
  end do
!************************************************************************** 
! Determinant calculation (ensure counter-clockwise orientation)
!**************************************************************************
  allocate(det(nel))
  do k=1,nel
    do i=2,3
       x(i)=z(1,mm(i,k))-z(1,mm(1,k))
       y(i)=z(2,mm(i,k))-z(2,mm(1,k))
    end do
    det(k)=x(2)*y(3)-y(2)*x(3)
    if (det(k)<0) then
      !call swap(mm(3,k),mm(2,k))
      temp = mm(2,k); mm(2,k) = mm(3,k); mm(3,k) = temp
    end if
  end do
!**************************************************************************
! nra: integer matrix (3,nel) where nra(i,k) is a reference number 
! associated to the i-th edge of the k-th element of the mesh

! nrv: integer matrix (3,nel) where nrv(i,k) is a reference number 
! associated to the i-th vertex of the k-th element of the mesh
!**************************************************************************
  allocate(nra(3,nel),nrv(3,nel),nrvg(nver),lfn(4))
  nra(1:3,1:nel)=0
  nrv(1:3,1:nel)=0
  nrvg(1:nver)=0  
  do j=ncarint+1,ndcar
    lfn(1:3)=0
    do i=1,2
      if (cver(i,j) == mm(1,celad(1,j))) then
        lfn(2)=1
      else if (cver(i,j) == mm(2,celad(1,j))) then
        lfn(3)=1
      else if (cver(i,j) == mm(3,celad(1,j))) then
        lfn(1)=1
      else if (cver(i,j) == mm(3,celad(1,j))) then
        lfn(2)=1
      end if      
      if (nrvg(cver(i,j)) == 0) then
        nrvg(cver(i,j))=bref(j)
      endif
    end do
    do k=1,3
      if (lfn(k) == 0) then
        nra(k,celad(1,j))=bref(j)
      end if  
    end do 
  end do   
  do k=1,nel
    nrv(1:3,k)=nrvg(mm(1:3,k))
  end do
!**************************************************************************
! nn: integer matrix (6,nel) where mm(i,k) is the global index of the i-th
! node of the k-th element of the mesh
!**************************************************************************
!  allocate(nn(6,nel),za(dim,3),ze(dim,3*nel))
!  nn(1:3,1:nel)=mm(1:3,1:nel)
!  na=0
!  do k=1,nel
!    za(1:3,1)=(z(1:3,mm(1,k))+z(1:3,mm(2,k)))/2
!    za(1:3,2)=(z(1:3,mm(2,k))+z(1:3,mm(3,k)))/2
!    za(1:3,3)=(z(1:3,mm(3,k))+z(1:3,mm(1,k)))/2
!bj: do j=1,3
!      do i=1,na
!        if (sum(abs(za(1:3,j)-ze(1:3,i)))<1e-6) then
!          nn(3+j,k)=nver+i
!          cycle bj
!        end if
!      end do
!    na=na+1
!    ze(1:3,na)=za(1:3,j)
!    nn(3+j,k)=nver+na
!    end do bj
!  end do
!**************************************************************************
! Write MUM file
!**************************************************************************
! Triangle data
!  if (output_mode == 1) then
!    nnoel=3
!    nnod=nver
!  else
!    nnoel=6
!    nnod=nver+na
!  endif

  nnoel=3
!  nnod=nver  
  npuel=3
  narel=3
  ncael=1

  lnv = npuel
  lne = narel
  lnf = ncael

lnn = lnv; nnod = nver !Assume Lagrange P1

print'(a,i9)','Global number of elements: ', nel
print'(a,i9)','Global number of nodes:    ', nnod
print'(a,i9)','Global number of vertices: ', nver
print'(a,i9)','Space dimension:           ', dim
print'(a,i9)','Local number of nodes:     ', lnn
print'(a,i9)','Local number of vertices:  ', lnv
print'(a,i9)','Local number of edges:     ', lne
print'(a,i9)','Local number of faces:     ', lnf
!print'(a)','EDGE REFERENCES ARE NOT READ!'
!print'(a)','ANSYS file loaded!'


  ! Print global data  
!  print*,'========================================'
!  print*,'        Global data for tetP1 mesh'
!  print*,'========================================'
!  print*,'Space dimension           : ',dim
!  print*,'Global number of vertices : ',nver
!  print*,'Global number of elements : ',nel
!  print*,'========================================'
  ! Write MUM for P1 elements
!  open(unit=15,file=file_out,form='unformatted')
!  rewind(15)
!  write(15) nel,nver,nver,dim,nnoel,npuel,narel,ncael
!  write(15) ((mm(i,k),i=1,3),k=1,nel), &
!     &      ((nra(i,k),i=1,3),k=1,nel), &
!     &      ((nrv(i,k),i=1,3),k=1,nel), &
!     &      ((z(i,j),i=1,2),j=1,nver)
!  write(15) ((nsd(i)),i=1,nel)
!  close(15)
!  print*,'========================================'
!  print*,'              Output files'
!  print*,'========================================'
!  print*,'File written: ' ,file_out
!  print*,'========================================'

end subroutine

!-----------------------------------------------------------------------
! load_ansys_tetra: loads an ANSYS tetrahedra mesh file
!-----------------------------------------------------------------------
subroutine load_ansys_tetra(iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, nn, mm, nrc, nra, nrv, z, nsd)
integer,                   intent(in)  :: iu
integer,                   intent(out) :: nel, nnod, nver, dim, lnn, lnv, lne, lnf
integer,      allocatable, intent(out) :: nn(:,:), mm(:,:), nrc(:,:), nra(:,:), nrv(:,:), nsd(:)
real(real64), allocatable, intent(out) :: z(:,:)

integer :: npuel, narel, ncael, st, ndcar, nref, nsref, fi, li, i, j, k, k_car, idx, &
bctype, ncarint, refname, temp
integer, allocatable :: cver(:,:), celad(:,:), contador(:), bref(:), nrvg(:), lfn(:)
real(real64) :: x(4), y(4), w(4)
real(real64), allocatable :: det(:)
character(20) :: aux(5)
character(60) :: idvar, var, hex_ver, hex_car, hex_elem, f_i, l_i, hextype

! Read data file '.msh' from ANSYS
! Read arguments
!  call read_arguments(file_in, file_out, output_mode)
! Open ansys file
!  open(unit=iu,file=file_in,form='formatted')
  rewind(iu)

!lee la linea que me indica si son tria o tetra
!defino lnn lnv....

!if (tri)
!read resto
!hago resto
!else tetra

! Comments
  do i=1,7
    read(iu,'(a)')
  end do
! Global dimensions
  read(iu,'(a)') var! idx,dim
  idx = 2
  dim = 3
! Hexadecimal reading
  read(iu,*) var,var,var,hex_ver,var
  read(iu,*) var,var,var,hex_car,var
  read(iu,*) var,var,var,hex_elem,var
! Global vertices
  read(hex_ver,'(z20)') nver
! Global faces
  read(hex_car,'(z20)') ndcar
! Global elements
  read(hex_elem,'(z20)') nel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!  write(len,*) 1+int(log(float(nel))/log(16.))
!  len='z'//trim(adjustl(len))  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
! Allocatables
  allocate(z(dim,nver),cver(dim,ndcar),celad(2,ndcar))
  allocate(bref(ndcar),nsd(nel))
! Initialization
  bref(1:ndcar)=0
  cver(1:dim,1:ndcar)=0
  celad(1:2,1:ndcar)=0
  nref=1
  nsref=1
  refname=1
!  print*,'========================================'
!  print*,'               References'
!  print*,'========================================'
! Loop to EOF (end of file)
bc: do while (.true.)
      read(iu,*,iostat=st) idvar
      if (is_iostat_end(st)) exit
      if (trim(idvar) /= '(45') then
         backspace(iu)
         read(iu,*) idvar,var,f_i,l_i
         read(f_i,'(z20)') fi
         read(l_i,'(z20)') li
         if (trim(idvar) == '(10') then
            ! Coordinates of vertices
            do j=fi,li
              read(iu,*) z(1:dim,j)
            end do
            read(iu,'(a)')          
         elseif (trim(idvar) == '(13') then
            backspace(iu)
            read(iu,*) idvar,var,f_i,l_i,hextype
            read(hextype,'(z20)') bctype
            ! Inside faces
            if (bctype == 2) then
                do j=fi,li
                   read(iu,*) aux(1:5)
                   read(aux(1),'(z20)') cver(1,j)
                   read(aux(2),'(z20)') cver(2,j)
                   read(aux(3),'(z20)') cver(3,j)
                   read(aux(4),'(z20)') celad(1,j)
                   read(aux(5),'(z20)') celad(2,j)
                   ncarint=li
                end do
            ! Frontier references                
            else
                do j=fi,li
                   read(iu,*) aux(1:5)
                   read(aux(1),'(z20)') cver(1,j)
                   read(aux(2),'(z20)') cver(2,j)
                   read(aux(3),'(z20)') cver(3,j)
                   read(aux(4),'(z20)') celad(1,j)
                   read(aux(5),'(z20)') celad(2,j)
                   bref(j)=nref
                end do
                nref=nref+1
            end if
            read(iu,'(a)')
            read(iu,'(a)')
         elseif (trim(idvar) == '(12') then
            ! Subdomains
            do j=fi,li
              nsd(j)=nsref
            end do
            nsref=nsref+1
         end if
      else
         backspace(iu)
         read(iu,*) idvar,var,f_i,l_i
         if (trim(f_i) /= 'interior' .and. trim(f_i) /= 'fluid') then
          ! Print references
!            print*,l_i(1:len_trim(l_i)-4),' Correspond to reference number ',refname
            print'(a)','Reference number '//trim(string(refname))//', associated with group named "'//l_i(1:len_trim(l_i)-4)//'"'
            refname=refname+1
         else  
            cycle bc
         end if 
      end if
    end do bc
  close(iu)
!  print*,'========================================'
!  print*,'Calculating ... '
!  print*,'========================================'  
! mm: integer matrix (4,nel) where mm(i,k) is the global index of the i-th
! node of the k-th element of the mesh
  allocate(mm(4,nel),contador(nel))
  mm(1:4,1:nel)=0
  contador(1:nel)=0
  do j=1,ncarint
    do i=1,2
      k_car=celad(i,j)
      if (contador(k_car)==0) then
        mm(1:3,k_car)=cver(1:3,j)        
        contador(k_car)=1
      else if (contador(k_car)==1) then
        do k=1,3
          if (cver(k,j) /= mm(1,k_car) .and. &
              cver(k,j) /= mm(2,k_car) .and. &
              cver(k,j) /= mm(3,k_car)) then
           mm(4,k_car) = cver(k,j)
           contador(k_car)=2
           exit
          end if
        end do
      end if
    end do      
  end do
  do j=ncarint+1,ndcar
    k_car=celad(1,j)
    if (contador(k_car)==0) then
      mm(1:3,k_car)=cver(1:3,j)        
      contador(k_car)=1
    else if (contador(k_car)==1) then
      do k=1,3
        if (cver(k,j) /= mm(1,k_car) .and. &
            cver(k,j) /= mm(2,k_car) .and. &
            cver(k,j) /= mm(3,k_car)) then
         mm(4,k_car) = cver(k,j)
         contador(k_car)=2
         exit
        end if
      end do
    end if    
  end do
! Determinant calculation (ensure counter-clockwise orientation)
  allocate(det(nel))
  do k=1,nel
    do i=2,4
       x(i) = z(1,mm(i,k))-z(1,mm(1,k))
       y(i)=z(2,mm(i,k))-z(2,mm(1,k))
       w(i)=z(3,mm(i,k))-z(3,mm(1,k))
    end do
    det(k)=x(2)*y(3)*w(4)+x(4)*y(2)*w(3)+x(3)*y(4)*w(2) &
        -x(4)*y(3)*w(2)-x(3)*y(2)*w(4)-x(2)*y(4)*w(3)
    if (det(k)<0) then
!      call swap(mm(2,k),mm(3,k))
     temp = mm(2,k); mm(2,k) = mm(3,k); mm(3,k) = temp
    end if    
  end do
! nrc: integer matrix (4,nel) where nrc(i,k) is a reference number 
! associated to the i-th face of the k-th element of the mesh

! nrv: integer matrix (4,nel) where nrv(i,k) is a reference number 
! associated to the i-th vertex of the k-th element of the mesh
  allocate(nrc(4,nel),nrv(4,nel),nrvg(nver),lfn(4))
  nrc(1:4,1:nel)=0
  nrv(1:4,1:nel)=0
  nrvg(1:nver)=0  
  do j=ncarint+1,ndcar
    lfn(1:4)=0
    do i=1,3
      if (cver(i,j) == mm(1,celad(1,j))) then
        lfn(4)=1
      else if (cver(i,j) == mm(2,celad(1,j))) then
        lfn(2)=1
      else if (cver(i,j) == mm(3,celad(1,j))) then
        lfn(3)=1
      else if (cver(i,j) == mm(4,celad(1,j))) then
        lfn(1)=1
      end if      
      if (nrvg(cver(i,j)) == 0) then
        nrvg(cver(i,j))=bref(j)
      endif
    end do
    do k=1,4
      if (lfn(k) == 0) then
        nrc(k,celad(1,j))=bref(j)
      end if  
    end do 
  end do   
  do k=1,nel
    nrv(1:4,k)=nrvg(mm(1:4,k))
  end do
  
!  if (output_mode /= 1) then
!**************************************************************************
! nn: integer matrix (10,nel) where mm(i,k) is the global index of the i-th
! node of the k-th element of the mesh
!**************************************************************************
!    allocate(nn(10,nel),za(dim,6),ze(dim,6*nel))
!    nn(1:4,1:nel)=mm(1:4,1:nel)
!    na=0
!    do k=1,nel
!      za(1:3,1)=(z(1:3,mm(1,k))+z(1:3,mm(2,k)))/2
!      za(1:3,2)=(z(1:3,mm(2,k))+z(1:3,mm(3,k)))/2
!      za(1:3,3)=(z(1:3,mm(3,k))+z(1:3,mm(1,k)))/2
!      za(1:3,4)=(z(1:3,mm(1,k))+z(1:3,mm(4,k)))/2
!      za(1:3,5)=(z(1:3,mm(2,k))+z(1:3,mm(4,k)))/2
!      za(1:3,6)=(z(1:3,mm(3,k))+z(1:3,mm(4,k)))/2
!  bj: do j=1,6
!        do i=1,na
!          if (sum(abs(za(1:3,j)-ze(1:3,i)))<1e-6) then
!            nn(4+j,k)=nver+i
!            cycle bj
!          end if
!        end do
!        na=na+1
!        ze(1:3,na)=za(1:3,j)
!        nn(4+j,k)=nver+na
!        end do bj
!      end do
!  endif
!**************************************************************************
! Write MUM file
!**************************************************************************
!  ! Tetrahedron data
!  if (output_mode == 1) then
!    nnoel=4
!    nnod=nver
!  else
!    nnoel=10
!    nnod=nver+na
!  endif
  npuel=4
  narel=6
  ncael=4

  ! Print global data  
!  print*,'========================================'
!  print*,'        Global data for mesh'
!  print*,'========================================'
!  print*,'Space dimension           : ',dim
!  print*,'Global number of vertices : ',nver
!  print*,'Global number of nodes    : ',nnod
!  print*,'Global number of elements : ',nel
!  print*,'========================================'

  ! Edge data (no se lee!, inutil para pb. electromagnéticos con fuente lineal)
  allocate(nra(6,nel))
  nra(1:6,1:nel) = 0

  ! Write MUM

!  if (output_mode == 1) then
!    write(15) nel,nver,nver,dim,nnoel,npuel,narel,ncael
!    write(15) ((mm(i,k),i=1,4),k=1,nel), &  
!       &      ((nrc(i,k),i=1,4),k=1,nel), &
!       &      ((nra(i,k),i=1,6),k=1,nel), &
!       &      ((nrv(i,k),i=1,4),k=1,nel), &
!       &      ((z(i,j),i=1,3),j=1,nver)
!    write(15) ((nsd(i)),i=1,nel)
!  elseif (output_mode == 2) then
!    write(15) nel,nnod,nver,dim,nnoel,npuel,narel,ncael
!    write(15) ((mm(i,k),i=1,4),k=1,nel), &
!       &      ((nn(i,k),i=1,10),k=1,nel), &    
!       &      ((nrc(i,k),i=1,4),k=1,nel), &
!       &      ((nra(i,k),i=1,6),k=1,nel), &
!       &      ((nrv(i,k),i=1,4),k=1,nel), &
!       &      ((z(i,j),i=1,3),j=1,nver)
!    write(15) ((nsd(i)),i=1,nel)
!  endif


!  if (output_mode /= 1) then
!    allocate(znod(dim,nnod))
!    znod(1:3,1:nver)=z(1:3,1:nver)
!    znod(1:3,nver+1:nnod)=ze(1:3,1:na)
!    zout=trim(file_out)//'_znod.dat'
!    open(unit=16,file=zout,form='unformatted')
!    rewind(16)
!    write(16) ((znod(i,j),i=1,3),j=1,nnod)
!    close(16)
!    print*,'File written : ',zout
!  endif

!  call reduce_bandwidth()

  lnv = npuel
  lne = narel
  lnf = ncael

lnn = lnv; nnod = nver !Assume Lagrange P1

print'(a,i9)','Global number of elements: ', nel
print'(a,i9)','Global number of nodes:    ', nnod
print'(a,i9)','Global number of vertices: ', nver
print'(a,i9)','Space dimension:           ', dim
print'(a,i9)','Local number of nodes:     ', lnn
print'(a,i9)','Local number of vertices:  ', lnv
print'(a,i9)','Local number of edges:     ', lne
print'(a,i9)','Local number of faces:     ', lnf

!print'(a)','EDGE REFERENCES ARE NOT READ!'
!print'(a)','ANSYS file loaded!'


!  if (output_mode == 1) then
!    print*,'========================================'  
!    print*,'Writing lagrange P1 tetrahedra mesh to MFM ... '
!    print*,'========================================'  
!    call savemum2mfm(file_out, 15, nel, nnod, nver, lnn, dim, lnv, lnn, lnf, &
!              nn, mm, nrc, nra, nrv, z, nsd)

!  elseif (output_mode == 2) then
!    print*,'========================================'  
!    print*,'Writing lagrange P2 tetrahedra mesh to MFM ... '
!    print*,'========================================'  
!    call savemum2mfm(file_out, 15, nel, nnod, nver, nnoel, dim, npuel, narel, ncael, &
!              nn, mm, nrc, nra, nrv, z, nsd)
!  elseif (output_mode == 3) then
!    print*,'========================================'  
!    print*,'Setting edge nodes ... '
!    print*,'========================================'  
!    call set_edge_nodes()
!    print*,'========================================'  
!    print*,'Writing Whitney P2 tetrahedra mesh to MFM ... '
!    print*,'========================================'  
!    call savelagr2edge(file_out, 15, nel, nnod, nver, lnn, dim, lnv, lnn, lnf, &
!              nn2, mm, nrc, nra, nrv, z, nsd)
!  elseif (output_mode == 4) then
!    print*,'========================================'  
!    print*,'Setting face nodes ... '
!    print*,'========================================'  
!    call set_face_nodes()
!    print*,'========================================'  
!    print*,'Writing Raviart-Thomas P2 tetrahedra mesh to MFM ... '
!    print*,'========================================'  
!    call savelagr2face(file_out, 15, nel, nnod, nver, lnn, dim, lnv, lnn, lnf, &
!              nn2, mm, nrc, nra, nrv, z, nsd)
!  endif

!  print*,'File written : ',file_out

!  ! Checking results
!  print*,'Checking nn'
!  do k=1,nel
!    print*,k,nn(1:10,k)
!  end do
!  print*,'Checking nrc'
!  do k=1,nel
!    print*,k,nrc(1:4,k)
!  end do
!  print*,'Checking nrv'
!  do k=1,nel
!    print*,k,nrv(1:4,k)
!  end do
!  print*,'Checking nsd'
!  do k=1,nel
!    print*,k,nsd(k)
!  end do
end subroutine

end module


