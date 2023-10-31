SUBROUTINE InitA(context,nprow,npcol,myrow,mycol,&
                &myrank, nprocs, total, atomfit, point, intervals,atomtype, numatom, &
                &atominv, ordatom, totmaxfit, invfit, kdim, xyz, dipxyz,fdenorm, &
                & maxfit,fitmaxidx,numfit)

!use mpi
implicit none

include 'mpif.h'

!input

!-------------------------------------------------
! point : numero di punto della griglia energetica
! intervals : array con gli intervalli energetici
!-------------------------------------------------

integer, intent(in) :: context,nprow,npcol,myrow,mycol
integer, intent(in) :: myrank, nprocs, totmaxfit
integer, intent(in) :: point, atomtype, numatom, atominv(numatom),total, &
                       & atomfit(numatom), ordatom(atomtype+1)
real*8, intent(in) :: intervals(point+1), fdenorm(total)
integer, intent(in) :: invfit(atomtype,totmaxfit)

integer, intent(in) :: maxfit(atomtype), fitmaxidx(numatom), numfit(atomtype)


integer, intent(inout) :: kdim(point+1), xyz, dipxyz(3)

! MPI
integer :: ierr, istat, fa, amode


character(len=4) :: nsym, occu, ieps, coef, DIPE
character(len=6) :: orbit, numbas, inpart
character(len=4), allocatable :: sym(:)

integer, allocatable :: numorb(:), virt1(:), nbasis(:), npart(:,:)
integer, allocatable :: nocc(:), nvir(:), symocc(:), symvir(:), kind(:)
integer, allocatable :: accop(:), indacc(:), symacc(:), dimeps(:)
integer, allocatable :: idxsymocc(:), idxsymvir(:)

real*8, allocatable :: energy(:), coefocc(:,:),coefvir(:,:)
real*8, allocatable :: froc(:), eps(:), coeftot(:)
real*8, allocatable :: vectx(:), vecty(:), vectz(:)
real*8, allocatable :: dipx(:), dipy(:), dipz(:), dip(:,:)
real*8, allocatable :: magx(:), magy(:), magz(:), mag(:,:)
real*8, allocatable :: dipvelx(:), dipvely(:), dipvelz(:), dipvel(:,:)
real*8, allocatable :: vectxa(:), vectya(:), vectza(:)
real*8, allocatable :: vectxb(:), vectyb(:), vectzb(:)
real*8, allocatable :: vectxb2(:), vectyb2(:), vectzb2(:)



integer :: numsym, occl, fd, db, fma, fdv
integer :: i, j, dim, p , l , m, q, k, z, dimecc
integer :: dimbas, dimocc, dimvir, tempi, ncount, maxbasis
integer :: II, III, IIII, JJ, JJJ, JJJJ, KK, INSIG
integer :: tempocc, tempvir, dimdip
integer :: status(MPI_STATUS_SIZE)
integer (KIND=MPI_OFFSET_KIND) :: offset, size

real*8 :: t1, t2, dt, tempr, dE
real*8, parameter :: eVha = 27.2113961

OPEN (UNIT=103, FILE= 'coef_eig', STATUS='OLD', ACTION='READ',IOSTAT=istat)
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

coefocc = 0.
coefvir = 0.

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

deallocate(sym)

CLOSE (103)

dipxyz = 0

if (myrank == 0) then
   READ(*,*) DIPE
   READ(*,*) xyz
   READ(*,*) dipxyz(1:xyz)
   READ(*,*) accop(1:numsym)
endif
CALL MPI_BCAST(accop,numsym,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(xyz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(dipxyz,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

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
            dE = eps(dimeps(symacc(m+indacc(i)))+j) -eps(dimeps(i)+occl)
            p = 0
            DO k = 1, point
               IF (intervals(k)/eVha <= dE .and. dE <intervals(k+1)/eVha) then
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
allocate (vectx(maxbasis))
allocate (vecty(maxbasis))
allocate (vectz(maxbasis))
allocate (vectxa(maxbasis))
allocate (vectya(maxbasis))
allocate (vectza(maxbasis))
allocate (vectxb(maxbasis))
allocate (vectyb(maxbasis))
allocate (vectzb(maxbasis))
allocate (vectxb2(maxbasis))
allocate (vectyb2(maxbasis))
allocate (vectzb2(maxbasis))
allocate (dipx(dimecc))
allocate (dipy(dimecc))
allocate (dipz(dimecc))
allocate (magx(dimecc))
allocate (magy(dimecc))
allocate (magz(dimecc))
allocate (dipvelx(dimecc))
allocate (dipvely(dimecc))
allocate (dipvelz(dimecc))


OPEN (UNIT=105, FILE= 'dipole', STATUS='OLD', ACTION='READ', IOSTAT=istat)
   READ (105,*) orbit, dimdip
   allocate (dip(dimdip,3))
   READ (105,*) dip
CLOSE (105)

OPEN (UNIT=110, FILE= 'magnetic', STATUS='OLD', ACTION='READ', IOSTAT=istat)
   READ (110,*) orbit, dimdip
   allocate (mag(dimdip,3))
   READ (110,*) mag
CLOSE (110)

OPEN (UNIT=111, FILE= 'dipvel', STATUS='OLD', ACTION='READ', IOSTAT=istat)
   READ (111,*) orbit, dimdip
   allocate (dipvel(dimdip,3))
   READ (111,*) dipvel
CLOSE (111)

z = 1
DO i = 1, numsym
   DO l = 1, virt1(i) - 1
      occl = virt1(i) - l
      q = 0
      tempocc = idxsymocc(i) + occl
      DO m = 1, accop(i)
         
         DO JJ = 1, nbasis(symacc(m+indacc(i)))
            vectx(JJ) = 0.0
            vecty(JJ) = 0.0
            vectz(JJ) = 0.0
            vectxa(JJ) = 0.0
            vectya(JJ) = 0.0
            vectza(JJ) = 0.0
            vectxb(JJ) = 0.0
            vectyb(JJ) = 0.0
            vectzb(JJ) = 0.0
            JJJ = npart(symacc(m+indacc(i)),JJ)
            DO II = 1, nbasis(i)
               III = npart(i, II)
               IIII = max(III,JJJ)
               JJJJ = min(III,JJJ)
               INSIG = SIGN( 1 , III - JJJ)
               KK = (IIII*(IIII-1))/2+JJJJ
               vectx(JJ) = vectx(JJ) + coefocc(II,tempocc) * dip(KK,1)
               vecty(JJ) = vecty(JJ) + coefocc(II,tempocc) * dip(KK,2)
               vectz(JJ) = vectz(JJ) + coefocc(II,tempocc) * dip(KK,3)
               vectxa(JJ) = vectxa(JJ) + INSIG * coefocc(II,tempocc) * mag(KK,1)
               vectya(JJ) = vectya(JJ) + INSIG * coefocc(II,tempocc) * mag(KK,2)
               vectza(JJ) = vectza(JJ) + INSIG * coefocc(II,tempocc) * mag(KK,3)
               vectxb(JJ) = vectxb(JJ) + INSIG * coefocc(II,tempocc) * dipvel(KK,1)
               vectyb(JJ) = vectyb(JJ) + INSIG * coefocc(II,tempocc) * dipvel(KK,2)
               vectzb(JJ) = vectzb(JJ) + INSIG * coefocc(II,tempocc) * dipvel(KK,3)
            ENDDO
         ENDDO

         DO j = virt1(symacc(m+indacc(i))), numorb(symacc(m+indacc(i)))
            dE = eps(dimeps(symacc(m+indacc(i)))+j) -eps(dimeps(i)+occl)
            p = 0
            DO k = 1, point
               IF (intervals(k)/eVha <= dE .and. dE <intervals(k+1)/eVha) then
                  nocc(z) = occl
                  nvir(z) = j
                  symvir(z) = symacc(m+indacc(i))
                  symocc(z) = i
                  energy(z) = dE
                  kind(z) = k
                  
                  dipx(z) = 0.0
                  dipy(z) = 0.0
                  dipz(z) = 0.0
                  magx(z) = 0.0
                  magy(z) = 0.0
                  magz(z) = 0.0
                  dipvelx(z) = 0.0
                  dipvely(z) = 0.0
                  dipvelz(z) = 0.0                  
                  tempvir = idxsymvir(symvir(z)) + j - virt1(symvir(z)) + 1

!         DO JJ = 1, nbasis(symacc(m+indacc(i)))
!            vectxb2(JJ) = 0.0
!            vectyb2(JJ) = 0.0
!            vectzb2(JJ) = 0.0
!            JJJ = npart(symacc(m+indacc(i)),JJ)
!            DO II = 1, nbasis(i)
!               III = npart(i, II)
!               IIII = max(III,JJJ)
!               JJJJ = min(III,JJJ)
!               KK = (IIII*(IIII-1))/2+JJJJ
!               vectxb2(JJ) = vectxb2(JJ) + coefvir(JJ,tempvir) * dipvel(KK,1)
!               vectyb2(JJ) = vectyb2(JJ) + coefvir(JJ,tempvir) * dipvel(KK,2)
!               vectzb2(JJ) = vectzb2(JJ) + coefvir(JJ,tempvir) * dipvel(KK,3)
!            ENDDO
!         ENDDO







                  DO JJ = 1, nbasis(symvir(z))

                     JJJ = npart(i,JJ)
!                     write(*,*) 'JJJ', JJJ

                     III = npart(symvir(z), JJ)
!                     write(*,*) 'III', III

                     IIII = max(III,JJJ)
                     JJJJ = min(III,JJJ)
                     KK = (IIII*(IIII-1))/2+JJJJ

                     dipx(z) = dipx(z) + coefvir(JJ,tempvir) * vectx(JJ)
                     dipy(z) = dipy(z) + coefvir(JJ,tempvir) * vecty(JJ)
                     dipz(z) = dipz(z) + coefvir(JJ,tempvir) * vectz(JJ)                     
                     magx(z) = magx(z) + coefvir(JJ,tempvir) * vectxa(JJ)
                     magy(z) = magy(z) + coefvir(JJ,tempvir) * vectya(JJ)
                     magz(z) = magz(z) + coefvir(JJ,tempvir) * vectza(JJ)
 

!dipvelx(z) = vectxb(JJ) + vectxb2(JJ)
!dipvely(z) = vectyb(JJ) + vectyb2(JJ)
!dipvelz(z) = vectzb(JJ) + vectzb2(JJ)

                     dipvelx(z) = dipvelx(z) + (coefvir(JJ,tempvir)*vectxb(JJ))
                     dipvely(z) = dipvely(z) + (coefvir(JJ,tempvir)*vectyb(JJ))
                     dipvelz(z) = dipvelz(z) + (coefvir(JJ,tempvir)*vectzb(JJ))



                  ENDDO

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
deallocate (vectx,vecty,vectz,dip,mag,dipvel)
deallocate (vectxa,vectya,vectza,vectxb,vectyb,vectzb)
deallocate (vectxb2,vectyb2,vectzb2)

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

         tempr = dipx(j)
         dipx(j) = dipx(j+1)
         dipx(j+1) = tempr

         tempr = dipy(j)
         dipy(j) = dipy(j+1)
         dipy(j+1) = tempr

         tempr = dipz(j)
         dipz(j) = dipz(j+1)
         dipz(j+1) = tempr

         tempr = magx(j)
         magx(j) = magx(j+1)
         magx(j+1) = tempr

         tempr = magy(j)
         magy(j) = magy(j+1)
         magy(j+1) = tempr

         tempr = magz(j)
         magz(j) = magz(j+1)
         magz(j+1) = tempr

         tempr = dipvelx(j)
         dipvelx(j) = dipvelx(j+1)
         dipvelx(j+1) = tempr

         tempr = dipvely(j)
         dipvely(j) = dipvely(j+1)
         dipvely(j+1) = tempr

         tempr = dipvelz(j)
         dipvelz(j) = dipvelz(j+1)
         dipvelz(j+1) = tempr

      endif
   enddo
enddo

OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)

if (myrank == 0) then
   WRITE (200,*) 'num ecc.', dimecc
   WRITE (200,*) '   K         Ek         SymI NumI  SymF NumF' 
   do i = 1, dimecc
      write(200,201) kind(i),energy(i),symocc(i),nocc(i),symvir(i),nvir(i) 
   enddo
endif

CLOSE (200)

offset = 0
CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)
call MPI_FILE_OPEN(MPI_COMM_WORLD,'dipole.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fd,ierr)
if (myrank == 0) then

   CALL MPI_FILE_WRITE_AT(fd, offset, dipx(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
   offset = dimecc*db
   CALL MPI_FILE_WRITE_AT(fd, offset, dipy(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
   offset = offset*2
   CALL MPI_FILE_WRITE_AT(fd, offset, dipz(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
endif
call MPI_FILE_CLOSE(fd,ierr)

deallocate (dipx, dipy, dipz)

offset = 0
CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)
call MPI_FILE_OPEN(MPI_COMM_WORLD,'magnetic.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fma,ierr)
if (myrank == 0) then
   CALL MPI_FILE_WRITE_AT(fma, offset, magx(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
   offset = dimecc*db
   CALL MPI_FILE_WRITE_AT(fma, offset, magy(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
   offset = offset*2
   CALL MPI_FILE_WRITE_AT(fma, offset, magz(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
endif
call MPI_FILE_CLOSE(fma,ierr)

deallocate (magx, magy, magz)

!write(*,*) 'dipvelx(1)', dipvelx(1), 'energy(1)', energy(1)

do i = 1, dimecc
   dipvelx(i) = dipvelx(i)/energy(i)
   dipvely(i) = dipvely(i)/energy(i)
   dipvelz(i) = dipvelz(i)/energy(i)
enddo   

offset = 0
CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)
call MPI_FILE_OPEN(MPI_COMM_WORLD,'dipvel.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fdv,ierr)
if (myrank == 0) then

   CALL MPI_FILE_WRITE_AT(fdv, offset, dipvelx(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
   offset = dimecc*db
   CALL MPI_FILE_WRITE_AT(fdv, offset, dipvely(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
   offset = offset*2
   CALL MPI_FILE_WRITE_AT(fdv, offset, dipvelz(1), dimecc,&
   & MPI_DOUBLE_PRECISION,status, ierr)
endif
call MPI_FILE_CLOSE(fdv,ierr)

deallocate (dipvelx, dipvely, dipvelz)



CALL CalcA(context,nprow,npcol,myrow,mycol,&
&myrank,nprocs,dimecc,atomtype,totmaxfit,numatom&
&,total,numsym,maxbasis,dimocc,dimvir,kind,atominv,ordatom,nbasis&
&,npart,invfit,atomfit,idxsymocc,idxsymvir,nocc,nvir,symocc&
&,symvir,virt1,coefocc,coefvir,fdenorm,maxfit,fitmaxidx,numfit)


kdim = 0
do i = 1, point
   kdim(i+1) = kdim(i)
   do j = 1, dimecc
      if (kind(j) == i) then
         kdim(i+1) = kdim(i+1)+1
      endif
   enddo
enddo

deallocate (virt1,nbasis,npart)
deallocate (nocc,nvir,symocc,symvir,energy,kind)
deallocate (coefocc, coefvir,idxsymocc,idxsymvir)

201 format (1x,i5,2x,E15.8,2x,i2,i5,4x,i2,i5)

END SUBROUTINE InitA
