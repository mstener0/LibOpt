Subroutine CalcMatrixAIF(myrank, nprocs, total, atomfit, point, intervals, atomtype, numatom, &
                       &atominv, ordatom, totmaxfit, invfit)

use mpi
implicit none

!input

!-------------------------------------------------
! point : numero di punto della griglia energetica
! intervals : array con gli intervalli energetici
!-------------------------------------------------

integer, intent(in) :: myrank, nprocs, totmaxfit
integer, intent(in) :: point, atomtype, numatom, atominv(numatom), total, &
                       & atomfit(numatom), ordatom(atomtype+1)
real*8, intent(in) :: intervals(point+1)
integer, intent(in) :: invfit(atomtype,totmaxfit)

! MPI
integer :: ierr, istat


character(len=4), parameter :: fita='fita', fitb='fitb'
character(len=4) :: nsym, occu, ieps, coef, DIP, fiti
character(len=8) :: atompair
character(len=2) :: la, lb
character(len=6) :: orbit, numbas, inpart, nbp, limit
character(len=4), allocatable :: sym(:)

integer, allocatable :: numorb(:), virt1(:), nbasis(:), npart(:,:), nbptr(:)
integer, allocatable :: nocc(:), nvir(:), symocc(:), symvir(:), kind(:)
integer, allocatable :: accop(:), indacc(:), symacc(:), dimeps(:), block(:)
integer, allocatable :: idxsymocc(:), idxsymvir(:), invpart(:,:)
integer, allocatable :: idxrig(:)

real*8, allocatable :: energy(:), MatrixAk(:,:), coefocc(:,:), coefvir(:,:)
real*8, allocatable :: froc(:), eps(:), coeftot(:)

integer :: ja, jb, basa, basb, idxfit, maxecc, Rig, Col
integer :: numsym, occl, idxbas(numatom), dimblo, tempp, tempq
integer :: i, j, dim, p , l , m, q, k, z, dimecc, N, in
integer :: dimbas, dimocc, dimvir, tempi, ncount, maxbasis
integer :: tempocc, tempvir

real*8 :: t1, t2, dt, tempr, dE, pairfit, den, rapr
real*8, parameter :: eVha = 27.2113961

OPEN (UNIT=103, FILE= 'coef_eig', STATUS='OLD', ACTION='READ', IOSTAT=istat)
READ (103,*) nsym, numsym 


allocate (sym(numsym))
allocate (numorb(numsym))
allocate (virt1(numsym))
allocate (nbasis(numsym))
allocate (accop(numsym))
allocate (indacc(numsym))
allocate (dimeps(numsym))


dim = 0
dimbas = 0
do i = 1, numsym
   READ (103,*) sym(i)
   READ (103,*) numbas, nbasis(i)
   READ (103,*) orbit, numorb(i)
   dim = dim + numorb(i)
   dimbas = dimbas + nbasis(i)
   allocate (froc(numorb(i)))
   READ (103,*) occu, froc
   virt1(i) = numorb(i) + 1
   do j = 1, numorb(i)
      if (froc(j) == 0.) then
         virt1(i) = j
         exit
      endif
   enddo
   deallocate (froc)
enddo

maxbasis = nbasis(1)
do i = 2, numsym
   IF (nbasis(i) > maxbasis) maxbasis = nbasis(i)
enddo

allocate (npart(numsym,maxbasis))
allocate (eps(dim))

npart = 0


p = 0
q = 0
dimeps = 0
do i = 1, numsym
   READ (103,*) sym(i)
   READ (103,*) ieps
   READ (103,*) eps(p+1:p+numorb(i))
   READ (103,*) inpart
   READ (103,*) npart(i,1:nbasis(i))
   dimeps(i) = p
   p = p + numorb(i)
   q = q + nbasis(i)
enddo

allocate (idxsymocc(numsym))
allocate (idxsymvir(numsym))

dimvir = 0
dimocc = 0
do i = 1, numsym
   idxsymocc(i) = dimocc
   idxsymvir(i) = dimvir
   dimocc = dimocc + (virt1(i) - 1)
   dimvir = dimvir + (numorb(i) - virt1(i) + 1)
enddo
 
allocate (coefocc(0:maxbasis,dimocc))
allocate (coefvir(0:maxbasis,dimvir))

!if (myrank==0) write(*,*) 'dimocc', dimocc
!if (myrank==0) write(*,*) 'dimvir', dimvir
!if (myrank==0) write(*,*) 'numorb', numorb
!if (myrank==0) write(*,*) 'virt1', virt1



coefocc = 0.
coefvir = 0.

!if (myrank==0) write(*,*) 'a'
!if (myrank==0) write(*,*) 'idxsymocc', idxsymocc
!if (myrank==0) write(*,*) 'idxsymvir', idxsymvir
!if (myrank==0) write(*,*) 'maxbasis', maxbasis

do i = 1, numsym
   allocate (coeftot(numorb(i)*nbasis(i)))

   READ (103,*) sym(i)
   READ (103,*) coef
 
   READ (103,*) coeftot

   p = 0
   do j = idxsymocc(i)+1, idxsymocc(i)+(virt1(i)-1)
      coefocc(1:(nbasis(i)),j) = coeftot(p+1:p+nbasis(i))
      p = p + nbasis(i)
   enddo
   do j = idxsymvir(i)+1, idxsymvir(i)+(numorb(i)-virt1(i)+1)
      coefvir(1:(nbasis(i)),j) = coeftot(p+1:p+nbasis(i))
      p = p + nbasis(i)
   enddo

   deallocate (coeftot)
enddo


if (myrank==0) write(*,*) 'idxsymocc', idxsymocc
if (myrank==0) write(*,*) 'idxsymvir', idxsymvir

deallocate(sym)


CLOSE (103)

if (myrank == 0) then
   READ(*,*) DIP
   READ(*,*) accop(1:numsym)
endif
CALL MPI_BCAST(accop,numsym,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

dim = 0
indacc = 0
do i = 1, numsym
   indacc(i) = dim
   dim = dim + accop(i)   
enddo

allocate (symacc(dim))

if (myrank == 0) then
   READ(*,*) symacc(1:dim)
endif
CALL MPI_BCAST(symacc,dim,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

dimecc = 0
DO i = 1, numsym
   DO l = 1, virt1(i) - 1
      occl = virt1(i) - l
      q = 0
      DO m = 1, accop(i)
         DO j = virt1(symacc(m+indacc(i))), numorb(symacc(m+indacc(i)))
            dE = eps(dimeps(symacc(m+indacc(i)))+j) - eps(dimeps(i)+occl)
            p = 0
            DO k = 1, point
               IF (intervals(k)/eVha <= dE .and. dE < intervals(k+1)/eVha) then
                  dimecc = dimecc + 1
                  q = 1
                  p = 1
                  EXIT
               ENDIF
            ENDDO
            IF ( p == 0 ) EXIT
         ENDDO
      ENDDO
      IF ( q == 0 ) EXIT
   ENDDO
ENDDO

if (myrank==0) then
   WRITE (*,*) 'ntotecc', dimecc
endif

allocate (nocc(dimecc))
allocate (nvir(dimecc))
allocate (symocc(dimecc))
allocate (symvir(dimecc))
allocate (energy(dimecc))
allocate (kind(dimecc))

z = 1
DO i = 1, numsym
   DO l = 1, virt1(i) - 1
      occl = virt1(i) - l
      q = 0
      DO m = 1, accop(i)
         DO j = virt1(symacc(m+indacc(i))), numorb(symacc(m+indacc(i)))
            dE = eps(dimeps(symacc(m+indacc(i)))+j) - eps(dimeps(i)+occl)
            p = 0
            DO k = 1, point
               IF (intervals(k)/eVha <= dE .and. dE < intervals(k+1)/eVha) then
                  nocc(z) = occl
                  nvir(z) = j
                  symvir(z) = symacc(m+indacc(i))
                  symocc(z) = i
                  energy(z) = dE
                  kind(z) = k
                  z = z + 1
                  q = 1
                  p = 1
                  EXIT
               ENDIF
            ENDDO
            IF ( p == 0 ) EXIT
         ENDDO
      ENDDO
      IF ( q == 0 ) EXIT
   ENDDO
ENDDO

deallocate (numorb,accop,indacc,dimeps,eps,symacc)

do i = 1, dimecc
   do j = 1, dimecc-1
      if (energy(j).GT.energy(j+1)) then

         tempr = energy(j)
         energy(j) = energy(j+1)
         energy(j+1) = tempr

         tempi = nocc(j)
         nocc(j) = nocc(j+1)
         nocc(j+1)=tempi

         tempi = nvir(j)
         nvir(j) = nvir(j+1)
         nvir(j+1)=tempi  

         tempi = symocc(j)
         symocc(j) = symocc(j+1)
         symocc(j+1)=tempi

         tempi = symvir(j)
         symvir(j) = symvir(j+1)
         symvir(j+1)=tempi

         tempi = kind(j)
         kind(j) = kind(j+1)
         kind(j+1)=tempi

      endif
   enddo
enddo

if (myrank==0) then
   READ(*,*) limit, maxecc
endif
CALL MPI_BCAST(maxecc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if (myrank==0) then
   write(*,*) 'limit', maxecc
endif

N = dimecc/maxecc
den = maxecc
rapr = dimecc/den

if (rapr > N) then
   N = N + 1
endif

allocate (block(N+1))
block(1) = 0
block(N+1) = dimecc

dimblo = dimecc/N

do i = 2, N
   block(i) = block(i-1) + dimblo
   in = block(i)
   do j = in, dimecc
      if (kind(j) /= kind(j+1)) EXIT
      block(i) = block(i) +1
   enddo
enddo
      
OPEN (UNIT=104, FILE= 'pairfit', STATUS='OLD', ACTION='READ', IOSTAT=istat)
allocate (nbptr(atomtype+1))
READ (104,*) nbp, nbptr

idxbas = 1
do i = 1, numatom-1
   idxbas(i+1) = nbptr(atominv(i)+1) - nbptr(atominv(i)) + idxbas(i)
enddo

dim = 0
do i = 1, atomtype
   dim = dim + ((nbptr(i+1)-nbptr(i))*(ordatom(i+1)-ordatom(i)))
enddo

allocate (invpart(numsym,dim))

invpart = 0
do i = 1, numsym
   do j = 1, nbasis(i)
      k = npart(i,j)
      invpart(i,k) = j     
   enddo
enddo   

deallocate(npart)

tempi = total/nprocs
allocate (idxrig(0:nprocs))
idxrig = 0

do i = 1, nprocs-1
   idxrig(i) = idxrig(i-1) + tempi
enddo
idxrig(nprocs) = total   

!write (*,*) 'idxrig', idxrig
tempi = idxrig(myrank+1) - idxrig(myrank)
!write(*,*) 'tempi', tempi

do i = 1, N

   allocate (MatrixAk(tempi ,block(i+1)-block(i)))
   MatrixAk = 0.

   ncount = 0
   readloop0: do
      READ(104,*) atompair, ja, jb, fiti
      if (ja==0 .and. jb==0) EXIT
      ncount = ncount + 1
      if (ja==jb) then

         readloop1: do
            READ(104,*) la, lb, basa, basb
            ncount = ncount + 1
            if (basa==0 .and. basb==0) EXIT
            basa = basa - nbptr(atominv(ja)) + idxbas(ja)
            basb = basb - nbptr(atominv(jb)) + idxbas(jb)

            readloop2: do
               READ(104,*) idxfit, pairfit
               ncount = ncount + 1
               if (idxfit==0 .and. pairfit==0.) EXIT

               j = invfit(atominv(ja),idxfit)
               if (j == 0) CYCLE

               Rig = j + atomfit(ja)
               if (Rig <= idxrig(myrank) .or. Rig > idxrig(myrank+1)) CYCLE
               Rig = Rig - idxrig(myrank)

               do l = block(i)+1, block(i+1)
                  Col = l - block(i)
                  tempocc = idxsymocc(symocc(l)) + nocc(l)
                  tempvir = idxsymvir(symvir(l)) + nvir(l) - virt1(symvir(l)) + 1

                  if (basa == basb) then
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basa)
                     if (k==0 .or. m==0) CYCLE
                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                  else
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basb)
                     if (k==0 .or. m==0) GOTO 10
                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))
       
10                   k = invpart(symvir(l),basb)
                     m = invpart(symocc(l),basa)
                     if (k==0 .or. m==0) CYCLE
                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                  endif           
               enddo

            enddo readloop2

         enddo readloop1
   
      else if (fiti .EQ. fita) then

         readloop1a: do
            READ(104,*) la, lb, basa, basb
            ncount = ncount + 1
            if (basa==0 .and. basb==0) EXIT
            basa = basa - nbptr(atominv(ja)) + idxbas(ja)
            basb = basb - nbptr(atominv(jb)) + idxbas(jb)

            readloop2a: do
               READ(104,*) idxfit, pairfit
               ncount = ncount + 1
               if (idxfit==0 .and. pairfit==0.) EXIT

               j = invfit(atominv(ja),idxfit)
               if (j == 0) cycle

               Rig = j  + atomfit(ja)
               if (Rig <= idxrig(myrank) .or. Rig > idxrig(myrank+1)) CYCLE
               Rig = Rig - idxrig(myrank)

               do l = block(i)+1, block(i+1)
                  Col = l - block(i)
                  tempocc = idxsymocc(symocc(l)) + nocc(l)
                  tempvir = idxsymvir(symvir(l)) + nvir(l) - virt1(symvir(l))+1

                  k = invpart(symvir(l),basa)
                  m = invpart(symocc(l),basb)
                  if (k==0 .or. m==0) GOTO 11
                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

11                k = invpart(symvir(l),basb)
                  m = invpart(symocc(l),basa)
                  if (k==0 .or. m==0) CYCLE
                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

               enddo

            enddo readloop2a
      
         enddo readloop1a
   
      else if (fiti .EQ. fitb) then

         readloop1b: do
            READ(104,*) la, lb, basa, basb
            ncount = ncount + 1
            if (basa==0 .and. basb==0) EXIT
            basa = basa - nbptr(atominv(ja)) + idxbas(ja)
            basb = basb - nbptr(atominv(jb)) + idxbas(jb)

            readloop2b: do
               READ(104,*) idxfit, pairfit
               ncount = ncount + 1
               if (idxfit==0 .and. pairfit==0.) EXIT

               j = invfit(atominv(jb),idxfit)
               if (j == 0) cycle

               Rig = j + atomfit(jb)
               if (Rig <= idxrig(myrank) .or. Rig > idxrig(myrank+1)) cycle
               Rig = Rig - idxrig(myrank)

               do l = block(i)+1, block(i+1)
                  Col = l - block(i)
                  tempocc = idxsymocc(symocc(l)) + nocc(l)
                  tempvir = idxsymvir(symvir(l)) + nvir(l) - virt1(symvir(l))+1

                  k = invpart(symvir(l),basa)
                  m = invpart(symocc(l),basb)
                  if (k==0 .or. m==0) GOTO 12
                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

12                k = invpart(symvir(l),basb)
                  m = invpart(symocc(l),basa)
                  if (k==0 .or. m==0) CYCLE
                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

               enddo

            enddo readloop2b
      
         enddo readloop1b

      endif

   enddo readloop0

   if (i == 1) then
   do in = 0, nprocs-1
      if (myrank == in) then
         write(*,*) 'myrank', myrank
         write(*,*) 'MatrixAk', MatrixAk(:,1)
      endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   enddo
   endif


   deallocate (MatrixAk)
   REWIND(104)
   BACKSPACE(104)



enddo

CLOSE(104)

deallocate (invpart,virt1,nbasis)
deallocate (nbptr,block)
deallocate (nocc,nvir,symocc,symvir,energy,kind)
deallocate (coefocc, coefvir,idxsymocc,idxsymvir)
deallocate (idxrig)

9002 format(2x,1e22.14,1x,1e22.14,1x,1e22.14,1x,1e22.14,1x,1e22.14)

End Subroutine CalcMatrixAIF

