SUBROUTINE CalcA(context,nprow,npcol,myrow,mycol,&
&myrank,nprocs,dimecc,atomtype,totmaxfit,numatom&
&,total,numsym,maxbasis,dimocc,dimvir,kind,atominv,ordatom,nbasis&
&,npart,invfit,atomfit,idxsymocc,idxsymvir,nocc,nvir,symocc&
&,symvir,virt1,coefocc,coefvir,fdenorm,maxfit,fitmaxidx,numfit)

implicit none

include 'mpif.h'

!input

integer, intent(in) :: context,nprow,npcol,myrow,mycol

integer, intent(in) :: myrank, nprocs, dimecc, atomtype, totmaxfit,&
                       & numatom, total, numsym, maxbasis, dimocc, dimvir

integer, intent(in) :: atominv(numatom), atomfit(numatom), ordatom(atomtype+1),&
                       & kind(dimecc), nbasis(numsym), idxsymocc(numsym), &
                       & idxsymvir(numsym), nocc(dimecc), nvir(dimecc), &
                       & symocc(dimecc), symvir(dimecc), virt1(numsym)  

integer, intent(in) :: invfit(atomtype,totmaxfit), npart(numsym,maxbasis)

real*8, intent(in) :: coefocc(0:maxbasis,dimocc), coefvir(0:maxbasis,dimvir),&
                      & fdenorm(total)

integer, intent(in) :: maxfit(atomtype),fitmaxidx(numatom),numfit(atomtype)

!variables

integer :: ierr, istat

character(len=4), parameter :: fita='fita', fitb='fitb'

character(len=4) :: fiti
character(len=8) :: atompair
character(len=2) :: la, lb
character(len=6) :: nbp, limit
character(len=5) :: fit3c
character(len=1) :: YN

integer, allocatable :: nbptr(:), block(:), idxrig(:)
integer, allocatable :: invpart(:,:)

real*8, allocatable :: MatrixAk(:,:)

integer :: ja, jb, basa, basb, idxfit, maxecc, Rig, Col
integer :: idxbas(numatom), dimblo
integer :: i, j, dim, l, m, k, N, in
integer :: tempi, ncount
integer :: tempocc, tempvir, fa

real*8 :: t1, t2, dt, pairfit, den, rapr

integer (KIND=MPI_OFFSET_KIND) :: size

YN = 'N'

if (myrank==0) then
   READ(*,*) limit, maxecc
   READ(*,*) fit3c, YN

endif
CALL MPI_BCAST(maxecc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(YN,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

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
      
size = 0

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixAk.out',MPI_MODE_RDONLY&
                  &,MPI_INFO_NULL,fa,ierr)

call MPI_FILE_GET_SIZE(fa, size, ierr)

call MPI_FILE_CLOSE(fa,ierr)

if (size == 0) then

OPEN (UNIT=104, FILE= 'pairfit', STATUS='OLD', ACTION='READ', IOSTAT=istat)
allocate (nbptr(atomtype+1))
READ (104,*) nbp, nbptr

OPEN (UNIT=109, FILE= 'pairfitextra', STATUS='OLD', ACTION='READ',IOSTAT=istat)
READ (109,*) nbp, nbptr

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

tempi = total/nprocs
allocate (idxrig(0:nprocs))

idxrig = 0

do i = 1, nprocs-1
   idxrig(i) = idxrig(i-1) + tempi
enddo
idxrig(nprocs) = total   

tempi = idxrig(myrank+1) - idxrig(myrank)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixAk.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fa,ierr)

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
                  tempvir = idxsymvir(symvir(l)) + nvir(l) - virt1(symvir(l)) +1

                  if (basa == basb) then
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basa)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))
                  else
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basb)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                     k = invpart(symvir(l),basb)
                     m = invpart(symocc(l),basa)

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

                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                  k = invpart(symvir(l),basb)
                  m = invpart(symocc(l),basa)

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

                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                  k = invpart(symvir(l),basb)
                  m = invpart(symocc(l),basa)

                  MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                    &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

               enddo

            enddo readloop2b
      
         enddo readloop1b

      endif

   enddo readloop0

   ncount = 0
   readloop0extra: do
      READ(109,*) atompair, ja, jb, fiti
      if (ja==0 .and. jb==0) EXIT
      ncount = ncount + 1
      if (fiti .EQ. fitb) then

         readloop1extraa: do
            READ(109,*) la, lb, basa, basb
            ncount = ncount + 1
            if (basa==0 .and. basb==0) EXIT
            basa = basa - nbptr(atominv(ja)) + idxbas(ja)
            basb = basb - nbptr(atominv(ja)) + idxbas(ja)

            readloop2extraa: do
               READ(109,*) idxfit, pairfit
               ncount = ncount + 1
               if (idxfit==0 .and. pairfit==0.) EXIT

               j = invfit(atominv(jb),idxfit)
               if (j == 0) CYCLE

               Rig = j + atomfit(jb)
               if (Rig <= idxrig(myrank) .or. Rig > idxrig(myrank+1)) CYCLE
               Rig = Rig - idxrig(myrank)

               do l = block(i)+1, block(i+1)
                  Col = l - block(i)
                  tempocc = idxsymocc(symocc(l)) + nocc(l)
                  tempvir = idxsymvir(symvir(l)) + nvir(l) - virt1(symvir(l)) +1

                  if (basa == basb) then
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basa)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))
                  else
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basb)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                     k = invpart(symvir(l),basb)
                     m = invpart(symocc(l),basa)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                  endif
               enddo

            enddo readloop2extraa

         enddo readloop1extraa

      else if (fiti .EQ. fita) then

         readloop1extrab: do
            READ(109,*) la, lb, basa, basb
            ncount = ncount + 1
            if (basa==0 .and. basb==0) EXIT
            basa = basa - nbptr(atominv(jb)) + idxbas(jb)
            basb = basb - nbptr(atominv(jb)) + idxbas(jb)

            readloop2extrab: do
               READ(109,*) idxfit, pairfit
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
                  tempvir = idxsymvir(symvir(l)) + nvir(l) - virt1(symvir(l)) +1

                  if (basa == basb) then
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basa)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))
                  else
                     k = invpart(symvir(l),basa)
                     m = invpart(symocc(l),basb)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                     k = invpart(symvir(l),basb)
                     m = invpart(symocc(l),basa)

                     MatrixAk(Rig,Col) = MatrixAk(Rig,Col) + &
                                       &(pairfit*coefvir(k,tempvir)*coefocc(m,tempocc))

                  endif
               enddo

            enddo readloop2extrab

         enddo readloop1extrab

      endif

   enddo readloop0extra

   do j = 1, tempi
      MatrixAk(j,:) = MatrixAk(j,:)*fdenorm(j+idxrig(myrank))
   enddo    

   CALL StoreA(total,tempi,N,myrank,nprocs,i,fa,block,idxrig,MatrixAk)

   deallocate (MatrixAk)
   REWIND(104)
   BACKSPACE(104)
   REWIND(109)
   BACKSPACE(109)

enddo

CALL MPI_FILE_CLOSE(fa,ierr)

CLOSE(104)
CLOSE(109)

deallocate(block,nbptr,idxrig)

IF (YN == 'Y') then
   call Calcpairfit3c(context,nprow,npcol,myrow,mycol,&
        &dimecc,atomtype,numatom&
        &,total,numsym,maxbasis,dimocc,dimvir,atominv&
        &,atomfit,idxsymocc,idxsymvir,nocc,nvir,symocc&
        &,symvir,virt1,coefocc,coefvir,fdenorm,maxfit,fitmaxidx,numfit,&
        & idxbas,dim,invpart)
ENDIF

deallocate(invpart)

endif !size

3001 format (1x,E15.8,2x,i3,1x,i3,1xE15.8,2x,i3,1x,i3,1xE15.8)

END SUBROUTINE CalcA
