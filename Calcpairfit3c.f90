Subroutine Calcpairfit3c(context,nprow,npcol,myrow,mycol,&
&dimecc,atomtype,numatom&
&,total,numsym,maxbasis,dimocc,dimvir,atominv&
&,atomfit,idxsymocc,idxsymvir,nocc,nvir,symocc&
&,symvir,virt1,coefocc,coefvir,fdenorm,maxfit,fitmaxidx,numfit,&
&idxbas,dim,invpart)

implicit none

include 'mpif.h'

!input

integer, intent(in) :: context,nprow,npcol,myrow,mycol,dim

integer, intent(in) :: dimecc, atomtype, &
                       & numatom, total, numsym, maxbasis, dimocc,dimvir

integer, intent(in) :: atominv(numatom), atomfit(numatom),&
                       & idxsymocc(numsym), &
                       & idxsymvir(numsym), nocc(dimecc), nvir(dimecc),&
                       & symocc(dimecc), symvir(dimecc), virt1(numsym)  

real*8, intent(in) :: coefocc(0:maxbasis,dimocc),coefvir(0:maxbasis,dimvir),&
                      & fdenorm(total)

integer, intent(in) :: maxfit(atomtype), fitmaxidx(numatom), numfit(atomtype),&
                       & idxbas(numatom), invpart(numsym,dim)

!MPI
external :: descinit, numroc, pdelset, INDXL2G
integer :: numroc, INDXL2G
integer :: ierr, info, fs3c, fak3c, fak
integer :: status(MPI_STATUS_SIZE)

!variables

real*8, parameter :: alpha= 1.0d0, beta= 0.0d0

integer :: istat

character(len=4), parameter :: fita='fita', fitb='fitb'

character(len=4) :: fiti
character(len=8) :: atompair
character(len=2) :: la, lb
character(len=6) :: nbp, limit

integer, allocatable :: nbptr(:), block(:), idxrig(:), numcyc(:)

real*8, allocatable :: S3C(:,:), S3Ca(:,:), S3Cb(:,:)
real*8, allocatable :: vecAB(:), pairfit3C(:,:), Ak3C(:,:), Ak(:,:)
real*8, allocatable :: coefmat(:,:)

integer, allocatable :: ibasa(:), ibasb(:)
integer :: cplbas


integer :: ja, jb, basa, basb, idxfit, maxecc, jatmp
integer :: dimbloa, cyc
integer :: i, j, l, m, k, N, in, idx, irig,im
integer :: tempi, ncount, db
integer :: tempocc, tempvir, fa
integer :: k2, m2, ilalb

real*8 :: coeff
real*8 :: t1, t2, dt


integer :: tempa, tempb

! Matrix
integer :: Mirig
integer :: Mloc, Nloc, Nloca, Nlocb, vecloc, cbloc,ncol, eccloc
integer :: descS3C(9), descS3Ca(9), descS3Cb(9), descvecAB(9)
integer :: descpf3C(9), descAk3C(9)
integer :: desccoef(9)

real*8 :: valamatfit
integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp


CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)

OPEN (UNIT=115, FILE= 'amatfit', STATUS='OLD', ACTION='READ',IOSTAT=istat)
allocate (nbptr(atomtype+1))
READ (115,*) nbp, nbptr


OPEN (UNIT=116, FILE= 'logout', STATUS='NEW', ACTION='WRITE',IOSTAT=istat)

Mloc = numroc( total, total/nprow, myrow, 0, nprow )

eccloc = numroc( dimecc, dimecc/npcol, mycol, 0, npcol )

call descinit(descAk3C, total, dimecc, total/nprow, dimecc/npcol,&
                & 0, 0, context, max(1,Mloc), info)

allocate ( Ak3C ( Mloc, eccloc ) )
Ak3C = 0.0d0



ncount = 0
jatmp = 0

readloop0: do   
   READ (115,*) atompair, ja, jb
   if (ja==0 .and. jb==0) EXIT

   if (myrow==0 .and. mycol==0) then
      write(116,*) 'ja', ja, 'jb', jb
   endif

   ncount = ncount + 1   


   cplbas = (nbptr(atominv(ja)+1)-nbptr(atominv(ja)))* &
          & (nbptr(atominv(jb)+1)-nbptr(atominv(jb)))

   irig = 0

   allocate ( ibasa (cplbas) )
   allocate ( ibasb (cplbas) )
   ibasa = 0
   ibasb = 0
    
   cbloc = numroc( cplbas, cplbas/npcol, mycol, 0, npcol )  

call descinit(descpf3C, total, cplbas, total/nprow, cplbas/npcol,&
                & 0, 0, context, max(1,Mloc), info)

allocate ( pairfit3C ( Mloc, cbloc ) )
pairfit3C = 0.0d0


   ncol = maxfit(atominv(ja)) + maxfit(atominv(jb))


   Nloc = numroc( ncol, ncol/npcol, mycol, 0, npcol )

   if (ja/=jatmp .and. ja /=2) then
   deallocate (S3Ca)
   endif


   if (ja/=jatmp) then
   Nloca = numroc( maxfit(atominv(ja)), maxfit(atominv(ja))/npcol, mycol, 0, npcol )
   call descinit(descS3Ca, total, maxfit(atominv(ja)), total/nprow, maxfit(atominv(ja))/npcol,&
                & 0, 0, context, max(1,Mloc), info)
   allocate ( S3Ca ( Mloc, Nloca ) )   
   S3Ca = 0.0d0
   CALL MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS3C.out',MPI_MODE_RDONLY&
                     &,MPI_INFO_NULL,fs3c,ierr)

   DO i = 1, Nloca
      col = INDXL2G(i, maxfit(atominv(ja))/npcol, mycol, 0, npcol)
!      col = INDXL2G(i, descS3Ca(5), mycol, 0, npcol)


      temp = (fitmaxidx(ja) + col - 1) * total
      IF (descS3Ca(9)>descS3Ca(5)) then
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp + rig -1) * db
         CALL MPI_FILE_READ_AT(fs3c, offset, S3Ca(1,i), descS3Ca(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
         offset = (temp + rig -1 +(descS3Ca(5)*nprow)) * db
         CALL MPI_FILE_READ_AT(fs3c,offset,S3Ca(1+descS3Ca(5),i),descS3Ca(9)-descS3Ca(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
      ELSE
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp+rig-1)*db
         CALL MPI_FILE_READ_AT(fs3c,offset,S3Ca(1,i),descS3Ca(5),&
                              & MPI_DOUBLE_PRECISION,status,ierr)
      ENDIF
   ENDDO

   call MPI_FILE_CLOSE(fs3c,ierr)   

   endif



   Nlocb = numroc( maxfit(atominv(jb)), maxfit(atominv(jb))/npcol, mycol, 0, npcol )


   vecloc = numroc( ncol, ncol/nprow, myrow, 0, nprow )



   call descinit(descS3C, total, ncol, total/nprow, ncol/npcol,&
                & 0, 0, context, max(1,Mloc), info)
!   call descinit(descS3Ca, total, maxfit(atominv(ja)), total/nprow, maxfit(atominv(ja))/npcol,&
!                & 0, 0, context, max(1,Mloc), info)
   call descinit(descS3Cb, total, maxfit(atominv(jb)), total/nprow, maxfit(atominv(jb))/npcol,&
                & 0, 0, context, max(1,Mloc), info)
   call descinit(descvecAB, ncol, 1, ncol/nprow, 1,&
                & 0, 0, context, max(1,vecloc), info)


   allocate ( S3C ( Mloc, Nloc ) )
!   allocate ( S3Ca ( Mloc, Nloca ) )
   allocate ( S3Cb ( Mloc, Nlocb ) )
   allocate ( vecAB ( vecloc ) )

S3C = 0.0d0
!S3Ca = 0.0d0
S3Cb = 0.0d0
vecAB = 0.0d0




   CALL MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS3C.out',MPI_MODE_RDONLY&
                     &,MPI_INFO_NULL,fs3c,ierr)

!   DO i = 1, Nloca
!      col = INDXL2G(i, maxfit(atominv(ja))/npcol, mycol, 0, npcol)
!      col = INDXL2G(i, descS3Ca(5), mycol, 0, npcol)


!      temp = (fitmaxidx(ja) + col - 1) * total
!      IF (descS3Ca(9)>descS3Ca(5)) then
!         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
!         offset = (temp + rig -1) * db
!         CALL MPI_FILE_READ_AT(fs3c, offset, S3Ca(1,i), descS3Ca(5),&
!                              & MPI_DOUBLE_PRECISION,status, ierr)
!         offset = (temp + rig -1 +(descS3Ca(5)*nprow)) * db
!         CALL MPI_FILE_READ_AT(fs3c,offset,S3Ca(1+descS3Ca(5),i),descS3Ca(9)-descS3Ca(5),&
!                              & MPI_DOUBLE_PRECISION,status, ierr)
!      ELSE
!         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
!         offset = (temp+rig-1)*db
!         CALL MPI_FILE_READ_AT(fs3c,offset,S3Ca(1,i),descS3Ca(5),&
!                              & MPI_DOUBLE_PRECISION,status,ierr)
!      ENDIF
!   ENDDO      

   DO i = 1, Nlocb

      col = INDXL2G(i, maxfit(atominv(jb))/npcol, mycol, 0, npcol)
!      col = INDXL2G(i, descS3Cb(5), mycol, 0, npcol)

      temp = (fitmaxidx(jb) + col - 1) * total
      IF (descS3Cb(9)>descS3Cb(5)) then
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp + rig -1) * db
         CALL MPI_FILE_READ_AT(fs3c, offset, S3Cb(1,i), descS3Cb(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
         offset = (temp + rig -1 +(descS3Cb(5)*nprow)) * db
         CALL MPI_FILE_READ_AT(fs3c,offset,S3Cb(1+descS3Cb(5),i),descS3Cb(9)-descS3Cb(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
      ELSE
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp+rig-1)*db
         CALL MPI_FILE_READ_AT(fs3c,offset,S3Cb(1,i),descS3Cb(5),&
                              & MPI_DOUBLE_PRECISION,status,ierr)
      ENDIF
   ENDDO

   call MPI_FILE_CLOSE(fs3c,ierr)

   call pdgeadd('N', total,maxfit(atominv(ja)),alpha,S3Ca,1,1,descS3Ca,beta,S3C,1,1,descS3C)
   call pdgeadd('N', total,maxfit(atominv(jb)),alpha,S3Cb,1,1,descS3Cb,beta,S3C,1,maxfit(atominv(ja))+1,descS3C)


!   deallocate (S3Ca, S3Cb)
   deallocate (S3Cb)  


   do i = 1, ncol
      do j = 1, numfit(atominv(ja))
         call pdelset(S3C, atomfit(ja)+j, i, descS3C, beta)
      enddo
      do j = 1, numfit(atominv(jb))
         call pdelset(S3C, atomfit(jb)+j, i, descS3C, beta)               
      enddo   
   enddo

   if (ja /=jatmp)then
   tempa = idxbas(ja) - nbptr(atominv(ja))
   endif

   jatmp = ja


   tempb = idxbas(jb) - nbptr(atominv(jb))

   readloop1: do
      READ(115,*) la, lb, basa, basb   
      ncount = ncount + 1
      if (basa==0 .and. basb==0) EXIT
      
      irig = irig + 1
      ibasa(irig) = basa + tempa
      ibasb(irig) = basb + tempb

      vecAB =0.0d0

      ilalb = 1

      readloop2: do 
         READ(115,*) idxfit, valamatfit
         ncount = ncount +1

         if (idxfit==0 .and. ilalb==1) then
            irig = irig - 1
            GOTO 100
         endif

         ilalb = ilalb + 1

         if (idxfit==0 .and. valamatfit==0.) EXIT
      
         call pdelset(vecAB, idxfit, 1, descvecAB, valamatfit)


      enddo readloop2

      CALL pdgemv('N',total,ncol,alpha,S3C,1,1,descS3C,vecAB,1,1,descvecAB,&
                    & 1,beta,pairfit3C,1,irig,descpf3C,1)

100   enddo readloop1

   deallocate (S3C, vecAB)

   if (irig ==0 ) then
      deallocate (ibasa, ibasb)
      deallocate (pairfit3C)
      CYCLE
   endif   

   if (irig < nprow) then
      if (myrow< irig) then
         Mirig = 1
      else
         Mirig = 0
      endif
   else
      Mirig = numroc( irig, irig/nprow, myrow, 0, nprow )
   endif   

   if (irig < nprow) then
      call descinit(desccoef, irig, dimecc, 1, dimecc/npcol,&
                   & 0, 0, context, max(1,Mirig), info)
   else
      call descinit(desccoef, irig, dimecc, irig/nprow, dimecc/npcol,&
                   & 0, 0, context, max(1,Mirig), info)
   endif

   allocate ( coefmat ( Mirig, eccloc ) )
   coefmat = 0.0d0   

   do j = 1, dimecc
      tempocc = idxsymocc(symocc(j)) + nocc(j)
      tempvir = idxsymvir(symvir(j)) + nvir(j) - virt1(symvir(j)) +1

      do i=1,irig

         k = invpart(symvir(j),ibasa(i))
         m = invpart(symocc(j),ibasb(i))

         k2 = invpart(symvir(j),ibasb(i))
         m2 = invpart(symocc(j),ibasa(i))

         call pdelset(coefmat,i,j,desccoef,coefvir(k,tempvir)*coefocc(m,tempocc)+&
                     &coefvir(k2,tempvir)*coefocc(m2,tempocc))

      enddo
    enddo   

   deallocate (ibasa, ibasb)

   CALL pdgemm('N','N',total,dimecc,irig,alpha,pairfit3C,1,1,descpf3C,coefmat,1,1,desccoef,&
                 & alpha,Ak3C,1,1,descAk3C)

   deallocate(coefmat)

   deallocate (pairfit3C)






!do j = 1, dimecc 
!      tempocc = idxsymocc(symocc(j)) + nocc(j)
!      tempvir = idxsymvir(symvir(j)) + nvir(j) - virt1(symvir(j)) +1
!
!    do i=1,irig
!
!      k = invpart(symvir(j),ibasa(i))
!      m = invpart(symocc(j),ibasb(i))
!
!      coeff = coefvir(k,tempvir)*coefocc(m,tempocc)
!
!      call pdgeadd('N',total,1,coeff,pairfit3C,1,i,descpf3C,alpha,Ak3C,1,j,descAk3C)
!
!      k = invpart(symvir(j),ibasb(i))
!      m = invpart(symocc(j),ibasa(i))
!
!      coeff = coefvir(k,tempvir)*coefocc(m,tempocc)
!
!      call pdgeadd('N',total,1,coeff,pairfit3C,1,i,descpf3C,alpha,Ak3C,1,j,descAk3C)
!
!   enddo
!enddo

!adeallocate (ibasa, ibasb)
!deallocate (nbptr)
!deallocate (pairfit3C)

enddo readloop0

   if (myrow==0 .and. mycol==0) then
      write(116,*) 'end'
   endif

CLOSE(116)



allocate ( Ak ( Mloc, eccloc ) )
Ak = 0.0d0

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixAk.out',MPI_MODE_RDWR&
                  &,MPI_INFO_NULL,fak,ierr)

DO i = 1, eccloc
   col = INDXL2G(i, dimecc/npcol, mycol, 0, npcol)

   temp = (col - 1) * total
   IF (descAk3C(9)>descAk3C(5)) then
      rig = INDXL2G(1, descAk3C(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_READ_AT(fak, offset, Ak(1,i), descAk3C(5),&
                           & MPI_DOUBLE_PRECISION,status, ierr)
      rig = INDXL2G(descAk3C(5)+1, descAk3C(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_READ_AT(fak,offset,Ak(1+descAk3C(5),i),descAk3C(9)-descAk3C(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      rig = INDXL2G(1, descAk3C(5), myrow, 0, nprow)
      offset = (temp+rig-1)*db
      CALL MPI_FILE_READ_AT(fak,offset,Ak(1,i),descAk3C(5),&
                            & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO


call pdgeadd('N', total,dimecc,alpha,Ak3C,1,1,descAk3C,alpha,Ak,1,1,descAk3C)

deallocate (Ak3C)

DO i = 1, eccloc
   col = INDXL2G(i, dimecc/npcol, mycol, 0, npcol)

   temp = (col - 1) * total
   IF (descAk3C(9)>descAk3C(5)) then
      rig = INDXL2G(1, descAk3C(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fak, offset, Ak(1,i), descAk3C(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)

      rig = INDXL2G(descAk3C(5)+1, descAk3C(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fak,offset,Ak(1+descAk3C(5),i),descAk3C(9)-descAk3C(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      rig = INDXL2G(1, descAk3C(5), myrow, 0, nprow)
      offset = (temp+rig-1)*db
      CALL MPI_FILE_WRITE_AT(fak,offset,Ak(1,i),descAk3C(5),&
                            & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO



call MPI_FILE_CLOSE(fak,ierr)

!call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixAk3C.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
!                  &,MPI_INFO_NULL,fak3c,ierr)
!
!DO i = 1, eccloc
!   col = INDXL2G(i, dimecc/npcol, mycol, 0, npcol)
!
!   temp = (col - 1) * total
!   IF (descAk3C(9)>descAk3C(5)) then
!      rig = INDXL2G(1, descAk3C(5), myrow, 0, nprow)
!      offset = (temp + rig -1) * db
!      CALL MPI_FILE_WRITE_AT(fak3c, offset, Ak3C(1,i), descAk3C(5),&
!                            & MPI_DOUBLE_PRECISION,status, ierr)
!
!      rig = INDXL2G(descAk3C(5)+1, descAk3C(5), myrow, 0, nprow)
!      offset = (temp + rig -1) * db
!      CALL MPI_FILE_WRITE_AT(fak3c,offset,Ak3C(1+descAk3C(5),i),descAk3C(9)-descAk3C(5),&
!                            & MPI_DOUBLE_PRECISION,status, ierr)
!   ELSE
!      rig = INDXL2G(1, descAk3C(5), myrow, 0, nprow)
!      offset = (temp+rig-1)*db
!      CALL MPI_FILE_WRITE_AT(fak3c,offset,Ak3C(1,i),descAk3C(5),&
!                            & MPI_DOUBLE_PRECISION,status,ierr)
!   ENDIF
!ENDDO
!
!call MPI_FILE_CLOSE(fak3c,ierr)



deallocate (Ak,nbptr)

end subroutine Calcpairfit3c

