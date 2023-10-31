Subroutine CalcMatrixL(context,nprow,npcol,myrow,mycol,total,numatom,atomtype,atominv,&
                      &atomfit,totmaxfit,invfit,fdenorm,lamb,totalmax,fitmaxidx)

!==========================================================!
! Purpose: Calc matrix L                                   !
! L = S^-1(F+Z)                                            !
!==========================================================!

use mpi

implicit none

!----------------------------------------------------------!
!input                                                     !
!                                                          !
!blacs variables: context, nprow, npcol, myrow, mycol      !
!total: number of selection fitting function               !
!numatom: number of total atoms                            !
!atomtype: number of type of atoms                         !
!totmaxfit: number of total fitting function               !
!atominv: array to pass from number of atom to type of atom!
!atomfit: array wirh first fitting function per atom       !
!invfit: array to select fitting function                  !
!fdenorm: array to denormalization of fitting function     !
!                                                          !
!----------------------------------------------------------!

integer, intent(in) :: context, nprow,npcol,myrow,mycol, total, numatom, atomtype, totmaxfit
integer, intent(in) :: totalmax, fitmaxidx(numatom)
integer, intent(in) :: atominv(numatom),atomfit(numatom),invfit(atomtype,totmaxfit)
real*8, intent(in) :: lamb
real*8, intent(inout) :: fdenorm(total)


! MPI
integer :: ierr, info, fh, fs, fs3c
integer :: status(MPI_STATUS_SIZE)
external :: descinit, numroc, pdelset, INDXL2G
integer :: numroc, INDXL2G

! Matrix
integer :: Nloc, Mloc, S3CNloc 
real*8, allocatable :: A(:,:), B(:,:)!, work(:)
real*8, allocatable :: S3C(:,:)
!integer, allocatable :: ipiv(:)
integer :: descMat(9), descS3C(9)
real*8 :: aval, val

character*9 :: atompair
character*4 :: fiti
integer :: i, j, k, l, p, q, ncount, db
integer :: ja, jb, nfaa, nfbb, istat
integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp

real*8 :: t1, t2, dt


! Computation of local matrix size
Mloc = numroc( total, total/nprow, myrow, 0, nprow )
Nloc = numroc( total, total/nprow, mycol, 0, npcol )

S3CNloc = numroc( totalmax, total/nprow, mycol, 0, npcol )

allocate( A( Mloc, Nloc ) )
allocate( B( Mloc, Nloc ) )

allocate( S3C(Mloc,S3CNloc))

! Descriptos
call descinit(descMat, total, total, total/nprow, total/nprow, 0, 0, context, max(1,Mloc), info)

call descinit(descS3C, total, totalmax, total/nprow, total/nprow, 0, 0, context, max(1,Mloc), info)

A = 0.0d0
B = 0.0d0 
fdenorm = 0.0d0

S3C = 0.0d0

! Read overlap fitting function matrix

OPEN (UNIT=100, FILE= 'sfit', STATUS='OLD', ACTION='READ', IOSTAT=istat)
do i = 1, numatom
   do j = 1, i
      read (100,*) atompair, ja, jb
      if ( ja == 0 .and. jb == 0 ) goto 1000
      read (100,*) fiti, nfaa, fiti, nfbb
      ncount = 0
      readloop0: do
         read (100, 6031) l, k, aval
         if ( k == 0 .or. l == 0 ) exit
         ncount = ncount + 1

         p = invfit(atominv(ja),k)
         q = invfit(atominv(jb),l)              

         if (ja/=jb) then
            if (p/=0) then
               call pdelset(S3C, p+atomfit(ja),l+fitmaxidx(jb), descS3C, aval) 
            endif
            if (q/=0) then
               call pdelset(S3C, q+atomfit(jb),k+fitmaxidx(ja), descS3C, aval)
            endif
         endif

         if (p==0 .or. q==0) CYCLE
         if (ja==jb .and. k == l) then
            fdenorm(p+atomfit(ja)) = SQRT(aval)
         endif
         call pdelset( A, p+atomfit(ja),q+atomfit(jb), descMat, aval)
         call pdelset( A, q+atomfit(jb),p+atomfit(ja), descMat, aval)
          
      enddo readloop0
   enddo
enddo


1000 CLOSE (100)

!Write overlap fitting function matrix in binary file

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fs,ierr)

CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)

DO i = 1, Nloc
   col = INDXL2G(i, total/nprow, mycol, 0, npcol)
   temp = (col - 1) * total
   IF (descMat(9)>descMat(5)) then
      rig = INDXL2G(1, descMat(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fs, offset, A(1,i), descMat(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)

      rig = INDXL2G(descMat(5)+1, descMat(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fs,offset,A(1+descMat(5),i),descMat(9)-descMat(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      rig = INDXL2G(1, descMat(5), myrow, 0, nprow)
      offset = (temp+rig-1)*db
      CALL MPI_FILE_WRITE_AT(fs,offset,A(1,i),descMat(5),&
                            & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO

call MPI_FILE_CLOSE(fs,ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS3C.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fs3c,ierr)

DO i = 1, S3CNloc
   col = INDXL2G(i, total/nprow, mycol, 0, npcol)
   temp = (col - 1) * total
   IF (descS3C(9)>descS3C(5)) then
      rig = INDXL2G(1, descS3C(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fs3c, offset, S3C(1,i), descS3C(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)

      rig = INDXL2G(descS3C(5)+1, descS3C(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fs3c,offset,S3C(1+descS3C(5),i),descS3C(9)-descS3C(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      rig = INDXL2G(1, descS3C(5), myrow, 0, nprow)
      offset = (temp+rig-1)*db
      CALL MPI_FILE_WRITE_AT(fs3c,offset,S3C(1,i),descS3C(5),&
                            & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO

call MPI_FILE_CLOSE(fs3c,ierr)

deallocate (S3C)

!Read matrix F

OPEN (UNIT=101, FILE= 'matrixF', STATUS='OLD', ACTION='READ', IOSTAT=istat)
do i = 1, numatom
   do j = 1, i
      read (101,*) atompair, ja, jb
      read (101,*) fiti, nfaa, fiti, nfbb
      ncount = 0
      readloop1: do
         read (101, 6031) l, k, aval
         if ( k == 0 .or. l == 0 ) exit
         ncount = ncount + 1

         p = invfit(atominv(ja),k)
         q = invfit(atominv(jb),l)
         if (p==0 .or. q==0) CYCLE
         call pdelset( B, p+atomfit(ja),q+atomfit(jb), descMat, aval)
         call pdelset( B, q+atomfit(jb),p+atomfit(ja), descMat, aval)

      enddo readloop1
   enddo
enddo

CLOSE (101)


!Read matrix Z and sum with matrix F

OPEN (UNIT=102, FILE= 'matrixZ', STATUS='OLD', ACTION='READ', IOSTAT=istat)
do i = 1, numatom
   do j = 1, i
      read (102,*) atompair, ja, jb
      read (102,*) fiti, nfaa, fiti, nfbb
      ncount = 0
      readloop2: do
         read (102, 6031) k, l, aval
         if ( k == 0 .or. l == 0 ) exit
         ncount = ncount + 1

         p = invfit(atominv(ja),k)
         q = invfit(atominv(jb),l)
         if (p==0 .or. q==0) CYCLE
         call pdelget('A', 'D', val, B, p+atomfit(ja),q+atomfit(jb), descMat)
         aval = aval + val
         call pdelset(B, p+atomfit(ja), q+atomfit(jb), descMat, aval*lamb)
         call pdelset(B, q+atomfit(jb), p+atomfit(ja), descMat, aval*lamb)

      enddo readloop2
   enddo
enddo
CLOSE (102)

! Linear system equations solver
!call pdgesv(total, total, A, 1, 1, descMat, ipiv, B, 1, 1, descMat, info )
call pdposv('U', total, total, A, 1, 1, descMat, B, 1, 1, descMat, info)

!Write matrix L in binary file

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixL.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fh,ierr)

CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)

DO i = 1, Nloc
   col = INDXL2G(i, total/nprow, mycol, 0, npcol)
   temp = (col - 1) * total
   IF (descMat(9)>descMat(5)) then
      rig = INDXL2G(1, descMat(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fh, offset, B(1,i), descMat(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)
      rig = INDXL2G(descMat(5)+1, descMat(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_WRITE_AT(fh,offset,B(1+descMat(5),i),descMat(9)-descMat(5),&
                            & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      rig = INDXL2G(1, descMat(5), myrow, 0, nprow)
      offset = (temp+rig-1)*db
      CALL MPI_FILE_WRITE_AT(fh,offset,B(1,i),descMat(5),&
                            & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO

call MPI_FILE_CLOSE(fh,ierr)

deallocate( A )
deallocate( B )

6031 format(1x,i5,1x,i5,2x,1e22.14)

End Subroutine CalcMatrixL

