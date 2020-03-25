module module_dataset_2414_fcnv
!-----------------------------------------------------------------------
!Universal Dataset Number 2414
!Name:   Analysis Data
!Status: Current
!Owner:  Simulation
!Revision Date: 3-OCT-1994
!-----------------------------------------------------------------------
!Record 1:        FORMAT(1I10)
!                 Field 1       -- Analysis dataset label
!Record 2:        FORMAT(40A2)
!                 Field 1       -- Analysis dataset name
!Record 3:        FORMAT (1I10)
!                 Field 1:      -- Dataset location
!                                   1:    Data at nodes
!                                   2:    Data on elements
!                                   3:    Data at nodes on elements
!                                   5:    Data at points
!Record 4:        FORMAT (40A2)
!                 Field 1:      -- ID line 1
!Record 5:        FORMAT (40A2)
!                 Field 1:      -- ID line 2
!
!
!Record 6:        FORMAT (40A2)
!                 Field 1:      -- ID line 3
!Record 7:        FORMAT (40A2)
!                 Field 1:      -- ID line 4
!Record 8:        FORMAT (40A2)
!                 Field 1:      -- ID line 5
!Record 9:        FORMAT (6I10)
!                 Field 1:      -- Model type
!                                   0:   Unknown
!                                   1:   Structural
!                                   2:   Heat transfer
!                                   3:   Fluid flow
!                 Field 2:      -- Analysis type
!                                   0:   Unknown
!                                   1:   Static
!                                   2:   Normal mode
!                                   3:   Complex eigenvalue first order
!                                   4:   Transient
!                                   5:   Frequency response
!                                   6:   Buckling
!                                   7:   Complex eigenvalue second order
!                                   9:   Static non-linear
!                                  10:   Craig-Bampton constraint modes
!                                  11:   Equivalent attachment modes
!                                  12:   Effective mass modes
!                                  13:   Effective mass matrix
!                                  14:   Effective mass matrix
!                 Field 3:      -- Data characteristic
!                                   0:   Unknown
!                                   1:   Scalar
!                                   2:   3 DOF global translation vector
!                                   3:   6 DOF global translation & rotation
!                                         vector
!                                   4:   Symmetric global tensor
!                                   6:   Stress resultants
!                 Field 4:      -- Result type
!                                   2:Stress
!                                   3:Strain
!                                   4:Element Force
!                                   5:Temperature
!                                   6:Heat Flux
!                                   7:Strain Energy
!                                   8:Displacement
!                                   9:Reaction Force
!                                  10:Kinetic Energy
!                                  11:Velocity
!                                  12:Acceleration
!                                  13:Strain Energy Density
!                                  14:Kinetic Energy Density
!                                  15:Hydrostatic Pressure
!                                  16:Heat Gradient
!                                  17:Code Check Value
!                                  18:Coefficient of Pressure
!                                  19:Ply Stress
!                                  20:Ply Strain
!                                  21:Failure Index for Ply
!                                  22:Failure Index for Bonding
!                                  23:Reaction Heat Flow
!                                  24:Stress Error Density
!                                  25:Stress Variation
!                                  27:Element Stress Resultant
!                                  28:Length
!                                  29:Area
!                                  30:Volume
!                                  31:Mass
!                                  32:Constraint Force
!                                  34:Plastic Strain
!                                  35:Creep Strain
!                                  36:Strain Energy Error Norm
!                                  37:Dynamic Stress At Nodes
!                                  38:Heat Transfer Coefficient
!                                  39:Temperature Gradient
!                                  40:Kinetic Energy Dissipation Rate
!                                  41:Strain Energy Error
!                                  42:Mass Flow
!                                  43:Mass Flux
!                                  44:Heat Flow
!                                  45:View Factor
!                                  46:Heat Load
!                                  47:Stress Component
!                                  48:Green Strain
!                                  49:Contact Forces
!                                  50:Contact Pressure
!                                  51:Contact Stress
!                                  52:Contact Friction Stress
!                                  53:Velocity Component
!                                  54:Heat Flux Component
!                                  55:Infrared Heat Flux
!                                  56:Diffuse Solar Heat Flux
!                                  57:Collimated Solar Heat Flux
!                                  58:Safety Factor
!                                  59:Fatigue Damage
!                                  60:Fatigue Damage With Direction
!                                  61:Fatigue Life
!                                  62:Quality Index
!                                  94:Unknown Scalar
!                                  95:Unknown 3DOF Vector
!                                  96:Unknown 6DOF Vector
!                                  97:Unknown Symmetric Tensor
!                                  98:Unknown General Tensor
!                                  99:Unknown Stress Resultant
!                                 101:Gap Thickness
!                                 102:Solid Layer (+ surface)
!                                 103:Solid Layer (- surface)
!                                 104:Total Solid Layer
!                                 105:Flow Vector at Fill
!                                 106:Bulk Flow Vector
!                                 107:Core Displacement
!                                 108:Layered Shear Strain Rate
!                                 109:Shear Stress
!                                 110:Heat Flux (+ surface)
!                                 111:Heat Flux (- surface)
!                                 112:Layered Temperature
!                                 113:Bulk Temperature
!                                 114:Peak Temperature
!                                 115:Temperature at Fill
!                                 116:Mass Density
!                                 117:Pressure
!                                 118:Volumetric Skrinkage
!                                 119:Filling Time
!                                 120:Ejection Time
!                                 121:No-flow Time
!                                 122:Weld Line Meeting Angle
!                                 123:Weld Line Underflow
!                                 124:Original Runner Diameter
!                                 125:Optimized Runner Diameter
!                                 126:Change in Runner Diameter
!                                 127:Averaged Layered Cure
!                                 128:Layered Cure
!                                 129:Cure Rate
!                                 130:Cure Time
!                                 131:Induction Time
!                                 132:Temperature at Cure
!                                 133:Percent Gelation
!                                 134:Part Heat Flux (+ surface)
!                                 135:Part Heat Flux (- surface)
!                                 136:Part-Wall Temperature (+ surface)
!                                 137:Part-Wall Temperature (- surface)
!                                 138:Part Ejection Time
!                                 139:Part Peak Temperature
!                                 140:Part Average Temperature
!                                 141:Parting Temperature (+ surface)
!                                 142:Parting Temperature (- surface)
!                                 143:Parting Heat Flux (- surface)
!                                 144:Parting Heat Flux (+ surface)
!                                 145:Wall Temperature Convergence
!                                 146:Wall Temperature (- surface)
!                                 147:Wall Temperature (+ surface)
!                                 148:Line Heat Flux
!                                 149:Line Pressure
!                                 150:Reynold's Number
!                                 151:Line Film Coefficient
!                                 152:Line Temperature
!                                 153:Line Bulk Temperature
!                                 154:Mold Temperature
!                                 155:Mold Heat Flux
!                                 156:Rod Heater Temperature
!                                 157:Rod Heater Flux
!                                 158:Original Line Diameter
!                                 159:Optimized Line Diameter
!                                 160:Change in Line Diameter
!                                 161:Air Traps
!                                 162:Weld Lines
!                                 163:Injection Growth
!                                 164:Temp Diff (Celcius)
!                                 165:Shear Rate
!                                 166:Viscosity
!                                 167:Percentage
!                                 168:Time
!                                 169:Flow Direction
!                                 170:Speed
!                                 171:Flow Rate
!                                 172:Thickness Ratio
!                                 301:Sound Pressure
!                                 302:Sound Power
!                                 303:Sound Intensity
!                                 304:Sound Energy
!                                 305:Sound Energy Density
!                                >1000:  User defined result type
!                 Field 5:      -- Data type
!                                   1:   Integer
!                                   2:   Single precision floating point
!                                   4:   Double precision floating point
!                                   5:   Single precision complex
!                                   6:   Double precision complex
!                 Field 6:      -- Number of data values for the data
!                                  component (NVALDC)
!Record 10:       FORMAT (8I10)
!                 Field 1:      -- Integer analysis type specific data (1-8)
!Record 11:       FORMAT (8I10)
!                 Field 1:      -- Integer analysis type specific data (9,10)
!Record 12:       FORMAT (6E13.5)
!                 Field 1:      -- Real analysis type specific data (1-6)
!Record 13:       FORMAT (6E13.5)
!                 Field 1:      -- Real analysis type specific data (7-12)
!Note: See chart below for specific analysis type information.
!Dataset class: Data at nodes
!Record 14:       FORMAT (I10)
!                 Field 1:      -- Node number
!Record 15:       FORMAT (6E13.5)
!                 Fields 1-N:   -- Data at this node (NDVAL real or complex values)
!                 Note: Records 14 and 15 are repeated for each node.
!Dataset class: Data at elements
!Record 14:       FORMAT (2I10)
!                 Field 1:      -- Element number
!                 Field 2:      -- Number Of data values For this element(NDVAL)
!Record 15:       FORMAT (6E13.5)
!                 Fields 1-N:   -- Data on element(NDVAL Real Or Complex Values)
!                 Note: Records 14 and 15 are repeated for all elements.
!Dataset class: Data at nodes on elements
!RECORD 14:       FORMAT (4I10)
!                 Field 1:      -- Element number
!                 Field 2:      -- Data expansion code (IEXP)
!                                  1: Data present for all nodes
!                                  2: Data present for only 1st node -All other
!                                     nodes the same.
!                 Field 3:      -- Number of nodes on elements (NLOCS)
!                 Field 4:      -- Number of data values per node (NVLOC)
!RECORD 15:       FORMAT (6E13.5)
!                 Fields 1-N:   -- Data Values At Node 1 (NVLOC Real Or
!                                  Complex Values)
!                 Note:  Records 14 And 15 Are repeated For each Element.
!                        For Iexp = 1 Record 15 Is repeated NLOCS Times
!                        For Iexp = 2 Record 15 appears once
!Dataset class: Data at points
!RECORD 14:       FORMAT (5I10)
!                 Field 1:      -- Element number
!                 Field 2:      -- Data expansion code (IEXP)
!                                  1: Data present for all points
!                                  2: Data present for only 1st point -All other
!                                     points the same.
!                 Field 3:      -- Number of points on elements (NLOCS)
!                 Field 4:      -- Number of data values per point (NVLOC)
!                 Field 5:      -- Element order
!RECORD 15:       FORMAT (6E13.5)
!                 Fields 1-N:   -- Data Values At point 1 (NVLOC Real Or
!                                  Complex Values)
!                 Note:  Records 14 And 15 Are repeated For each Element.
!                        For Iexp = 1 Record 15 Is repeated NLOC Times
!                        For Iexp = 2 Record 15 appears once
!-----------------------------------------------------------------------
use basicmod
use module_dataset_fcnv
use module_pmh_fcnv, only: pmh_mesh, field
use module_fe_database_pmh_fcnv, only: FEDB
implicit none

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! read: read dataset 2411
!-----------------------------------------------------------------------
subroutine read_2414(iu, pmh, npc, nfield, els_loc, dataset, infield, ca_opt, padval,  param)
 integer,               intent(in) :: iu    !unit number for unvfile
 type(pmh_mesh),     intent(inout) :: pmh   !PMH mesh
 integer,               intent(in) :: npc   !Piece number
 integer,            intent(inout) :: nfield   !Piece number
 integer, allocatable,  intent(in) :: els_loc(:,:) !Elements location ([elgroup, pos], element label)
 integer,               intent(in) :: dataset
 character(len=*),allocatable, intent(in) :: infield(:) ! Input field names
 logical,               intent(in) :: ca_opt
 real(real64),          intent(in) :: padval
 real(real64), optional,intent(in) :: param
 integer, dimension(:)             :: r9(6),r10(8)
 real(real64), dimension(:)        :: r12(6)
 integer                           :: ios, counter, i, prev_nel
 integer                           :: nparam
 real(real64)                      :: paramval
 integer :: lbl, dloc, ncomp, fidx, n_nod, n_el, iexp, n_nod_el,nval
 character(len=maxpath) :: name, aux
 real(real64), allocatable :: val(:), tempval(:), mat3(:,:,:)
 type(field),  allocatable :: auxfi(:)
 logical :: found = .true.

  ! Single param
  nparam = 1
  paramval = 0._real64
  counter = 0

  if(dataset == 2414) then
    read (unit=iu, fmt='(1I10)', iostat = ios) lbl   ! Analysis dataset label. Record1
      if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
    read (unit=iu, fmt='(1A80)', iostat = ios) name  ! Analysis dataset name. Record 2
      if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
    read (unit=iu, fmt='(1I10)', iostat = ios) dloc  ! Dataset location, 1:nodes, 2:elements. Record3
      if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
    nfield = lbl
  elseif(dataset == 55) then
    nfield = nfield + 1
    lbl = nfield
    write(name,fmt='(I10.10)') nfield
    name = 'field'//trim(name)
    dloc = 1 ! Data at nodes
  elseif(dataset == 57 .or. dataset == 56) then
    nfield = nfield + 1
    lbl = nfield
    write(name,fmt='(I10.10)') nfield
    name = 'field'//trim(name)
    dloc = 2 ! Data at elements
  else
    call error('dataset_2414/read, # Dataset '//trim(string(dataset))//' not allowed')
  endif

  ! Read id lines
  read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read') !Record4
  read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read') !Record5
  read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read') !Record6
  read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read') !Record7
  read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read') !Record8
  ! Record 9
  read (unit=iu, fmt='(6I10)', iostat = ios) r9
    if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))

  if(r9(2) == 3 .or. r9(2) == 3 .or. r9(2) == 7) call error('dataset_2414/read, # Multivaluated parameter not supported') ! Multiple values for param

  if(r9(3) == 0) then
    ncomp = 0
  elseif(r9(3) == 1) then
    ncomp = 1 ! Scalar
  elseif(r9(3) == 2) then
    ncomp = 3 ! 3 Comp. vector
  elseif(r9(3) == 3 .or. r9(3) == 4) then
    ncomp = 6 ! 6 Comp. simmetric global tensor
  elseif(r9(3) == 5) then
    ncomp = 9 ! 9 Comp. general tensor
  else
    call error('dataset_2414/read, # Wrong data characteristic')
  endif

  ! Change fieldname with de Code Aster option
  if(ca_opt) then ! Option Code aster
    if(is_ca_field_type(string(r9))) then
      name = get_ca_field_name(string(r9)) !name = get_ca_field_name(string(r9), name)
    endif
  endif


  if(r9(5) == 5 .or. r9(5) == 6) call error('dataset_2414/read, # Complex data not supported') ! Compex data

  if(dataset == 2414) then
  ! Integer analysis specific data. Record10 & record11
    read (unit=iu, fmt='(8I10)', iostat = ios) r10; if (ios /= 0) call error('dataset_2414/read')
    read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read')
  ! Real analysis specific data. Record12 & record13
    read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read')
    read (unit=iu, fmt=*, iostat = ios) aux; if (ios /= 0) call error('dataset_2414/read')
  else
  ! Integer analysis specific data. Record7 & record8
    read (unit=iu, fmt='(8I10)', iostat = ios) r10; if (ios /= 0) call error('dataset_2414/read')
    read (unit=iu, fmt='(6E13.5)', iostat = ios) r12; if (ios /= 0) call error('dataset_2414/read')
    nparam = r10(4)
    paramval = r12(1)
  endif

  if(allocated(val)) deallocate(val); allocate(val(ncomp))
  if(allocated(tempval)) deallocate(tempval); allocate(tempval(ncomp))
  val = 0._real64
  tempval = 0._real64

  ! If allocated input field names
  if(allocated(infield)) then
    found = .false.
    ! Check if field is in file
    do i=1, size(infield,1)
      if(trim(adjustl(infield(i))) == trim(adjustl(name))) found = .true.
    enddo
    ! If field not in field names list, skip.
    if(.not. found) then
      do
        if (is_dataset_delimiter(iu, back=.false.)) return
        read (unit=iu, fmt=*, iostat = ios) aux
        if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
      enddo
    endif
  endif

  if(dloc == 1) then ! Data at nodes
    ! PMH field structure allocation
    if(.not. allocated(pmh%pc(npc)%fi)) then
      fidx = 1
      allocate(pmh%pc(npc)%fi(fidx))
      call info('Reading node field "'//trim(adjustl(name))//'" with record6/9: '//trim(string(r9)))
    else
      found = .false.
      do i=1, size(pmh%pc(npc)%fi,1)
        if(trim(adjustl(pmh%pc(npc)%fi(i)%name)) == trim(adjustl(name))) then
          fidx = i
          found = .true.
        endif
      enddo
      if(.not. found) then
        nparam = 1
        fidx = size(pmh%pc(npc)%fi,1)+1
        if(allocated(auxfi)) deallocate(auxfi)
        allocate(auxfi(fidx))
        auxfi(1:size(pmh%pc(npc)%fi,1)) = pmh%pc(npc)%fi(:)
        call move_alloc(from=auxfi, to=pmh%pc(npc)%fi)
        call info('Reading node field "'//trim(adjustl(name))//'" with record6/9: '//trim(string(r9)))
      else
        call info('  Reading shot "'//trim(string(nparam))//'" with parameter: '//trim(string(paramval)))
      endif
    endif
    ! PMH Field structure: name, param and values initilization
    pmh%pc(npc)%fi(fidx)%name = trim(adjustl(name))

    if(present(param)) then
      call set(pmh%pc(npc)%fi(fidx)%param,param,nparam,fit=.true.)
    else
      call set(pmh%pc(npc)%fi(fidx)%param,paramval,nparam,fit=.true.)
    endif
    call alloc(mat3, ncomp, pmh%pc(npc)%nnod, 1)
    mat3 = reshape([(padval,i=1,pmh%pc(npc)%nnod*ncomp)], [ncomp, pmh%pc(npc)%nnod, 1])
    call set(pmh%pc(npc)%fi(fidx)%val, mat3, [(i,i=1,ncomp)], [(i,i=1,pmh%pc(npc)%nnod)], [nparam], fit=[.true.,.true.,.true.])
    do
      if (is_dataset_delimiter(iu, back=.true.)) exit
    ! Node or element number. Record14
      read (unit=iu, fmt='(1I10)', iostat = ios) n_nod
      if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
    ! Data. Record15
      read (unit=iu, fmt=*, iostat = ios) pmh%pc(npc)%fi(fidx)%val(:,n_nod,nparam)
      if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
      counter  = counter + 1
    end do

    if (pmh%pc(npc)%nnod /= counter) call warning('dataset_2414/read, Maximum node numbering, '//trim(string(pmh%pc(npc)%nnod))//&
    ', differs from number of field entries, '//trim(string(counter))//'; possible non-consecutive node numbering. Continue '//&
    'at your own risk!')    
    !)Number of nodes # Wrong number of values')

  elseif(dloc == 2 .or. dloc == 3) then ! Data at elements
    fidx = 0
    do
      prev_nel = 0
      if (is_dataset_delimiter(iu, back=.true.)) exit
      ! Node or element number. Record14
      read (unit=iu, fmt='(4I10)', iostat = ios) n_el, iexp, n_nod_el, nval

      if(nval>6) then !if(nval/=0) then
        n_nod_el = n_nod_el*(nval/6) ! Number of rows with values
      endif

      if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
      ! Count elements in previous groups
      do i=1, els_loc(1,n_el)-1
        if(FEDB(pmh%pc(npc)%el(i)%type)%tdim>0) prev_nel = prev_nel + pmh%pc(npc)%el(i)%nel
      enddo

      associate(elg => pmh%pc(npc)%el(els_loc(1,n_el)))
        ! Check if field is a serie
        found = .false.
        if(allocated(elg%fi)) then
          do i=1, size(elg%fi,1)
            if(trim(adjustl(elg%fi(i)%name)) == trim(adjustl(name))) then
              fidx = i
              found = .true.
            endif
          enddo
        endif

        if(.not. found .and. .not. allocated(elg%fi)) then
          fidx = 1
          nparam = 1
          allocate(elg%fi(fidx))
          ! PMH Field structure: name, param and values initilization
          elg%fi(fidx)%name = trim(adjustl(name))
          call info('Reading cell field "'//trim(adjustl(name))//'" with record6/9: '//trim(string(r9)))
          if(.not. allocated(elg%fi(fidx)%param)) allocate(elg%fi(fidx)%param(nparam))
          if(.not. allocated(elg%fi(fidx)%val)) allocate(elg%fi(fidx)%val(ncomp,elg%nel,nparam))
          elg%fi(fidx)%param(nparam) = paramval
          elg%fi(fidx)%val(:,:,nparam) = padval
        elseif(.not. found .and. (fidx == 0 .or. fidx > size(elg%fi,1))) then
          fidx = size(elg%fi,1)+1
          nparam = 1
          if(allocated(auxfi)) deallocate(auxfi)
          call move_alloc(from=elg%fi, to=auxfi)
          allocate(elg%fi(fidx))
          elg%fi(1:fidx-1) = auxfi(:)
          if(allocated(auxfi)) deallocate(auxfi)
          ! PMH Field structure: name, param and values initilization
          elg%fi(fidx)%name = trim(adjustl(name))
          call info('Reading cell field "'//trim(adjustl(name))//'" with record6/9: '//trim(string(r9)))
          if(.not. allocated(elg%fi(fidx)%param)) allocate(elg%fi(fidx)%param(nparam))
          if(.not. allocated(elg%fi(fidx)%val)) allocate(elg%fi(fidx)%val(ncomp,elg%nel,nparam))
          elg%fi(fidx)%param(nparam) = paramval
          elg%fi(fidx)%val(:,:,nparam) = padval
        endif

        if(nparam > size(elg%fi(fidx)%param,1) .or. nparam > size(elg%fi(fidx)%val,3)) then
          if(found) call info('  Reading shot "'//trim(string(nparam))//'" with parameter: '//trim(string(paramval)))
          if(present(param)) then
            call set(elg%fi(fidx)%param,param,nparam,fit=.true.)
          else
            call set(elg%fi(fidx)%param,paramval,nparam,fit=.true.)
          endif
          call alloc(mat3, ncomp, elg%nel, 1)
          mat3 = reshape([(padval, i=1,elg%nel*ncomp)], [ncomp, elg%nel, 1])
          call set(elg%fi(fidx)%val, mat3, [(i,i=1,ncomp)], [(i,i=1,elg%nel)], [nparam], fit=[.true.,.true.,.true.])
        endif

      ! Data. Record15
        if((dataset == 57 .or. dataset == 56) .and. iexp == 2) then ! 2: Data present for only 1st node
          read (unit=iu, fmt=*, iostat = ios) elg%fi(fidx)%val(:,els_loc(2,n_el),nparam)
          if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
        elseif(((dataset == 57 .or. dataset == 56) .and. iexp == 1) .or. dloc == 3) then ! 1: Data present for all nodes
          val = 0._real64
          do i=1,n_nod_el
            read (unit=iu, fmt='('//trim(string(ncomp))//'E13.5)', iostat = ios) tempval
            if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
            val = val + tempval
          enddo
          elg%fi(fidx)%val(:,els_loc(2,n_el),nparam) = val/n_nod_el
        else
          read (unit=iu, fmt=*, iostat = ios) elg%fi(fidx)%val(:,els_loc(2,n_el),nparam)
          if (ios /= 0) call error('dataset_2414/read, #'//trim(string(ios)))
        endif
      end associate
    end do

  else
    call error('dataset_2414/read, # Data at elements not implemented yet')
  endif

end subroutine

!-----------------------------------------------------------------------
! is_ca_field_type(): returns true if input string is a field type identifier
! Field identifiers (record6 (55-56-57) or record9 (2414) are a list of 6I10.
!-----------------------------------------------------------------------
function is_ca_field_type(str) result(res)
  character(len=*), intent(in) :: str
  logical                      :: res
  integer                      :: ios, ints(6)

  res = .false.
  if( trim(adjustl(lcase(str))) == 'depl' .or. &
    & trim(adjustl(lcase(str))) == 'vite' .or. &
    & trim(adjustl(lcase(str))) == 'acce' .or. &
    & trim(adjustl(lcase(str))) == 'temp' .or. &
    & trim(adjustl(lcase(str))) == 'vari_elno' .or. &
    & trim(adjustl(lcase(str))) == 'epsa_elno' .or. &
    & trim(adjustl(lcase(str))) == 'sief_elno' .or. &
    & trim(adjustl(lcase(str))) == 'pres' ) then
    res = .true.
  else
    if(word_count(str) == 6) then
      read(unit=str, fmt='(6I10)', iostat = ios) ints
      if(ios == 0) res = .true.
    endif
  endif

end function

!-----------------------------------------------------------------------
! get_ca_field_type(): returns Code Aster field name if input string is a field type identifier
! Field identifiers (record6 (55-56-57) or record9 (2414) are a list of 6I10.
!-----------------------------------------------------------------------
function get_ca_field_name(str, inname, verbose) result(fname)
  character(len=*),          intent(in) :: str
  character(len=*),optional, intent(in) :: inname
  logical,         optional, intent(in) :: verbose
  character(len=maxpath)                :: fname
  logical                               :: vbs
  integer                               :: ios, ints(6)

  vbs = .false.
  if(present(verbose)) vbs=verbose

  if(.true.) then ! Code aster flag
    select case(trim(adjustl(lcase(str))))
      case('depl') ! Displacement Code Aster field
        fname = 'DEPL'; return
      case('vite')      ! Velocity Code Aster field
        fname = 'VITE'; return
      case('acce')      ! Acceleration Code Aster field
        fname = 'ACCE'; return
      case('temp')      ! Temperature Code Aster field
        fname = 'TEMP'; return
      case('vari_elno') ! VARI_ELNO Code Aster field
        fname = 'VARI_ELNO'; return
      case('epse_elno') ! EPSA_ELNO Code Aster field
        fname = 'EPSA_ELNO'; return
      case('sief_elno') ! SIEF_ELNO Code Aster field
        fname = 'SIEF_ELNO'; return
      case('pres')      ! Pressure Code Aster field
        fname = 'PRES'; return
      case DEFAULT
        if(word_count(str) == 6) then
          read(unit=str, fmt=*, iostat = ios) ints
          if(ios == 0) then
            if(    all(ints == (/1,4,3, 8,2,6/))) then ! Displacement Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'DEPL')
              fname = 'DEPL'; return
            elseif(all(ints == (/1,4,3,11,2,6/))) then ! Velocity Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'VITE')
              fname = 'VITE'; return
            elseif(all(ints == (/1,4,3,12,2,6/))) then ! Acceleration Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'ACCE')
              fname = 'ACCE'; return
            elseif(all(ints == (/2,4,1, 5,2,1/))) then ! Temperature Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'TEMP')
              fname = 'TEMP'; return
            elseif(all(ints == (/1,4,3, 0,2,6/))) then ! VARI_ELNO Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'VARI_ELNO')
              fname = 'VARI_ELNO'; return
            elseif(all(ints == (/1,4,4, 3,2,6/))) then ! EPSA_ELNO Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'EPSA_ELNO')
              fname = 'EPSA_ELNO'; return
            elseif(all(ints == (/1,4,4, 2,2,6/))) then ! SIEF_ELNO Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'SIEF_ELNO')
              fname = 'SIEF_ELNO'; return
            elseif(all(ints == (/1,4,1,15,2,1/))) then ! Pressure Code Aster field
              if(vbs .and. present(inname)) &
                & call info('  Field name '//trim(inname)//' changed to '//'PRES')
              fname = 'PRES'; return
            else
              if(present(inname)) then
                fname = inname; return
              else
                fname = trim(adjustl(str)); return;
              endif
            end if
          endif
        endif
    end select
  end if

end function

!-----------------------------------------------------------------------
! get_ca_field_record9(): returns field type identifies if input string is a code aster field name
! Field identifiers (record6 (55-56-57) or record9 (2414) are a list of 6I10.
!-----------------------------------------------------------------------
function get_ca_field_record9(str, inrecord) result(ints)
  character(len=*),          intent(in) :: str
  integer,optional, intent(in)          :: inrecord(6)
  integer                               :: ints(6)
  integer                               :: ios

  if(.true.) then ! Code aster flag
    select case(trim(adjustl(lcase(str))))
      case('depl') ! Displacement Code Aster field
        ints = (/1,4,3, 8,2,6/); return
      case('vite')      ! Velocity Code Aster field
        ints = (/1,4,3,11,2,6/); return
      case('acce')      ! Acceleration Code Aster field
        ints = (/1,4,3,12,2,6/); return
      case('temp')      ! Temperature Code Aster field
        ints = (/2,4,1, 5,2,1/); return
      case('vari_elno') ! VARI_ELNO Code Aster field
        ints = (/1,4,3, 0,2,6/); return
      case('epsa_elno') ! EPSA_ELNO Code Aster field
        ints = (/1,4,4, 3,2,6/); return
      case('sief_elno') ! SIEF_ELNO Code Aster field
        ints = (/1,4,4, 2,2,6/); return
      case('pres')      ! Pressure Code Aster field
        ints = (/1,4,1,15,2,1/); return
      case DEFAULT
        if(word_count(str) == 6) then
          read(unit=str, fmt=*, iostat = ios) ints
          if(ios == 0) then
            if( all(ints == (/1,4,3, 8,2,6/)) .or. &  ! Displacement Code Aster field
              & all(ints == (/1,4,3,11,2,6/)) .or. &  ! Velocity Code Aster field
              & all(ints == (/1,4,3,12,2,6/)) .or. &  ! Acceleration Code Aster field
              & all(ints == (/2,4,1, 5,2,1/)) .or. &  ! Temperature Code Aster field
              & all(ints == (/1,4,3, 0,2,6/)) .or. &  ! VARI_ELNO Code Aster field
              & all(ints == (/1,4,4, 3,2,6/)) .or. &  ! EPSA_ELNO Code Aster field
              & all(ints == (/1,4,4, 2,2,6/)) .or. &  ! SIEF_ELNO Code Aster field
              & all(ints == (/1,4,1,15,2,1/))) then   ! Pressure Code Aster field
              return
            else
              if(present(inrecord)) then
                ints = inrecord; return
              else
                ints = (/0,0,0,0,0,0/)
              endif
            end if
          endif
        endif
    end select
  end if

end function


end module
