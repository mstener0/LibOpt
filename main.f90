 program main

!=======================!
! LibOpt    ALGORITHM   !
!=======================!


use mpi

implicit none


! MPI
integer :: ierr, nprocs, myrank, istat

! Benchmarks
real*8 :: t1, t2, dt 


!-----------------------------------------------------------------------------------!
!                                                                                   !
! atomtype: number type of atoms                                                    !
! numatom: number of atoms                                                          !
! totfit: number of fitting function for type of atoms                              !
! total: total number of fitting functions                                          !
! ordatom(atomtype+1): array to order number of atom per type                       !
! numfit(atomtype): array with number of fitting function per atom type             !
! indfit(totfit): array with index of selected fitting function                     !
! atomfit(numatom): array with first fitting function per atom                      !
! atominv(numatom): array to pass from number of atom to type of atom               !
! invfit: array to select fitting function                                          !
!                                                                                   !
!-----------------------------------------------------------------------------------!


integer :: atomtype, numatom, totfit, total, totmaxfit, numpoint, totalmax
integer, allocatable :: ordatom(:), numfit(:), indfit(:), atomfit(:), &
                      & atominv(:), maxfit(:), kdim(:), fitidx(:), bindex(:)&
                      &, fitmaxidx(:)

character*8 :: ATOM, NUM, ANALYSIS
character*3, allocatable :: atoms(:)
character*3 :: FIT, var, sim
character*6 :: KCYCLE,ENERGY,LAMBDA,DIPVEL
character*4 :: VINC, ITER
character*1 :: YN, IT, DV, ANYN

integer :: l, ind, i, j, point, k, xyz, dipxyz(3), numstep
integer :: temp
real*8 :: En, range, lamb, Ein, Efin, wi, delta
real*8, allocatable :: intervals(:), fdenorm(:), enpoint(:)

integer, allocatable :: invfit(:,:)

! BLACS/SCALAPACK
integer :: context, nprow, npcol, myrow, mycol, info, ndims, dims(2)
parameter ( ndims=2 )
external :: blacs_exit, blacs_gridexit, blacs_gridinfo, blacs_get
external :: blacs_gridinit, descinit, numroc, pdelset, pdlaprnt, indxl2g
integer :: numroc, indxl2g
real*8, parameter :: eVha = 27.2113961


! Initialize MPI environment
call MPI_Init(ierr)
call MPI_Comm_Size(mpi_comm_world,nprocs,ierr)
call MPI_Comm_Rank(mpi_comm_world,myrank,ierr)

! Initialize a default BLACS context and the processes grid
dims = 0
call MPI_Dims_Create( nprocs, ndims, dims, ierr)
nprow = dims(1)
npcol = dims(2)
call blacs_get( -1, 0, context )
call blacs_gridinit( context, 'Row-major', nprow, npcol)
call blacs_gridinfo( context, nprow, npcol, myrow, mycol )

!=============!
!READ INPUT   !
!=============!

t1 = MPI_wtime()

ATOM = 'ATOMTYPE'
NUM = 'NUMATOM'

if (myrank == 0) then
   READ(*,*) ATOM, atomtype
   READ(*,*) NUM, numatom
endif

CALL MPI_BCAST (atomtype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (numatom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

allocate (ordatom(atomtype+1))
allocate (atoms(atomtype))
allocate (numfit(atomtype))
allocate (maxfit(atomtype))

ind = 0
ordatom = 0

if (myrank == 0) then
   do i = 1, atomtype
      READ(*,*) atoms(i), l
      ind = ind + l
      ordatom(i+1) = ind
   enddo
endif
CALL MPI_BCAST (ordatom(:), atomtype+1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

if (myrank == 0) then
   do i = 1, atomtype
      WRITE(*,*) atoms(i), ordatom(i)+1, ordatom(i+1)
   enddo
endif
CALL MPI_BCAST (atoms(:), atomtype*3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

FIT = 'FIT'

if (myrank == 0) then
   READ(*,*) var
   if (var == FIT) then
      do i = 1, atomtype
         READ(*,*) sim
         if (sim == atoms(i)) then
            READ(*,*) maxfit(i),numfit(i)
         endif
      enddo
   endif
endif
CALL MPI_BCAST (numfit(:), atomtype, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (maxfit(:), atomtype, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

if (myrank == 0) then
   write(*,*) 'numfit', numfit
endif

totmaxfit = maxfit(1)
do i = 2, atomtype
   IF (maxfit(i) > totmaxfit) totmaxfit = maxfit(i)
enddo

totfit = 0
do i = 1, atomtype
   totfit = totfit + numfit(i)
enddo

allocate (indfit(totfit))
allocate (invfit(atomtype,totmaxfit))

if (myrank == 0) then
   READ(*,*) indfit(1:totfit)
endif
CALL MPI_BCAST (indfit, totfit, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

invfit = 0

do i = 1, numfit(1)
   k = indfit(i)
   invfit(1,k) = i
enddo

allocate (fitidx(atomtype))
fitidx = 0

if (atomtype > 1) then
   temp = numfit(1)
   do i = 2, atomtype
      fitidx(i) = temp
      do j = temp+1, temp + numfit(i)
         k = indfit(j)
         invfit(i,k) = j - temp
      enddo
      temp = temp + numfit(i) 
   enddo
endif

allocate (atomfit(numatom))
allocate (atominv(numatom))
allocate (fitmaxidx(numatom))
atomfit = 0
total = 0
totalmax = 0
fitmaxidx = 0

do i = 1, atomtype
   do j = ordatom(i)+1, ordatom(i+1)
      total = total + numfit(i)
      totalmax = totalmax + maxfit(i)
      atominv(j) = i
      if (j == numatom) exit
      atomfit(j+1) = numfit(i) + atomfit(j)
      fitmaxidx(j+1) = maxfit(i) + fitmaxidx(j)
   enddo
enddo   

if (myrank ==0) then
   write(*,*) 'indfit', indfit
endif

allocate (bindex (total))

if (myrank == 0) OPEN (UNIT=200, FILE= 'analysis', STATUS='NEW', ACTION='WRITE',IOSTAT=istat)

do i = 1, atomtype
   do j = ordatom(i)+1, ordatom(i+1)
      bindex(atomfit(j)+1:atomfit(j)+numfit(i)) = fitmaxidx(j) + &
      &  indfit(fitidx(i)+1:fitidx(i)+numfit(i))
   enddo
enddo

if (myrank ==0 ) then
   write(200,*) 'b index', total
   write(200,*) bindex
endif

if (myrank == 0) CLOSE (200)

deallocate (fitidx,bindex)!,fitmaxidx)



call MPI_Barrier( mpi_comm_world, ierr)

if (myrank == 0) then
   read(*,*) LAMBDA, lamb
endif

CALL MPI_BCAST (lamb, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

if (lamb<0. .or. lamb>1.) STOP 'LAMBDA <0 or >1 '

allocate (fdenorm(total))

if (myrank ==0 ) then
   write(*,*) 'total', total
   write(*,*) 'numatom', numatom
   write(*,*) 'atomtype', atomtype
   write(*,*) 'atominv', atominv
   write(*,*) 'atomfit', atomfit
   write(*,*) 'totmaxfit', totmaxfit
   do i = 1, atomtype 
   write(*,*) 'invfit', invfit(i,:)
   enddo   
endif

CALL CalcMatrixL(context,nprow,npcol,myrow,mycol,total,numatom,&
                &atomtype,atominv,atomfit,totmaxfit,invfit,fdenorm,lamb,&
                &totalmax,fitmaxidx)

! lettura e creazione della griglia energetica

if (myrank == 0) then
   read(*,*) KCYCLE
   read(*,*) En, point
endif

CALL MPI_BCAST (En, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (point, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


range = En/point

allocate (intervals(point+1))
intervals = 0.

do i = 1, point
   intervals(i+1) = Intervals(i) + range
enddo

if (myrank==0) then
   write (*,*) 'En', En
   write (*,*) 'point', point
endif

allocate(kdim(point+1))
kdim = 0

if (myrank ==0 ) then
   write(*,*) 'ordatom', ordatom
endif


CALL InitA(context,nprow,npcol,myrow,mycol,myrank,nprocs,total,atomfit,point,&
          &intervals,atomtype,&
          &numatom,atominv,ordatom,totmaxfit,invfit, kdim, xyz, &
          &dipxyz,fdenorm,maxfit,fitmaxidx,numfit)

deallocate (fitmaxidx)

call MPI_Barrier( mpi_comm_world, ierr)

intervals = intervals / eVha

if (mycol == 0 .and. myrow == 0 ) then
   READ(*,*) VINC, YN
endif
CALL MPI_BCAST (YN, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

if (mycol == 0 .and. myrow == 0 ) then
   READ(*,*) ENERGY
   READ(*,*) Ein, Efin, numstep
   READ(*,*) wi
   READ(*,*) delta
endif

CALL MPI_BCAST (Ein, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (Efin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (numstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (wi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (delta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

delta = delta/eVha

if (mycol == 0 .and. myrow == 0 ) then
   READ(*,*) DIPVEL, DV
endif
CALL MPI_BCAST (DV, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)


if (mycol == 0 .and. myrow == 0 ) then
   READ(*,*) ANALYSIS, ANYN
endif
CALL MPI_BCAST (ANYN, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)


!if (mycol == 0 .and. myrow == 0 ) then
!   READ(*,*) ANALYSIS
!   READ(*,*) numpoint
!endif
!
!CALL MPI_BCAST (numpoint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!
!allocate (enpoint(numpoint))
!
!if (mycol == 0 .and. myrow == 0 ) then
!   READ(*,*) enpoint(1:numpoint)
!endif
!
!enpoint = enpoint/eVha
!
!CALL MPI_BCAST (enpoint, numpoint, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)



if (mycol == 0 .and. myrow == 0 ) then
   READ(*,*) ITER, IT
endif
CALL MPI_BCAST (IT, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

t2 = MPI_wtime()

dt = t2 - t1
call MPI_Barrier( mpi_comm_world, ierr)
write(*,*) 'tempo integrali', dt

t1 = MPI_wtime()



if (IT == 'Y') then
CALL Calcalpit(context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,&
          &dipxyz,atomtype,numatom,atominv,atomfit,totmaxfit,   &
          &invfit,maxfit,intervals,Ein,Efin,numstep,wi,YN,ANYN)
!numpoint,&
!          & enpoint)

else
CALL Calcalp(context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,&
          &dipxyz,atomtype,numatom,atominv,atomfit,totmaxfit,   &
          &invfit,maxfit,intervals,Ein,Efin,numstep,wi,YN,ANYN,delta,DV)
!numpoint,&
!          & enpoint)
endif

t2 = MPI_wtime()

dt = t2 - t1
call MPI_Barrier( mpi_comm_world, ierr)
write(*,*) 'tempo tddft', dt


deallocate( kdim )
deallocate( indfit )
deallocate( atoms )
deallocate( numfit )
deallocate( ordatom )
deallocate( atomfit )
deallocate( fdenorm )
deallocate( maxfit )
deallocate( invfit )
deallocate( intervals )
!deallocate( enpoint )

! Close BLACS environment
call blacs_gridexit( context )
call blacs_exit( 1 )


! Close MPI environment if blacs_exit paramater is not equal zero
call MPI_Finalize(ierr)

6031 format(1x,i5,1x,i5,2x,1e22.14)

end program main

