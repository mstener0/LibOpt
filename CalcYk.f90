SUBROUTINE CalcYk(context,nprow,npcol,myrow,mycol,total,point,kdim)

use mpi

implicit none

!input
integer, intent(in) :: context,nprow,npcol,myrow,mycol, total, point
integer, intent(in) :: kdim(point+1)

integer, parameter :: double = selected_real_kind(14)
real(double), parameter :: eVha = 27.2113961
integer :: i, j, db

! MPI
integer :: ierr, fa, fs, fy
integer :: status(MPI_STATUS_SIZE)
! BLACS/SCALAPACK
integer :: info
external :: blacs_exit, blacs_gridexit, blacs_gridinfo, blacs_get, &
            &blacs_gridinit, blacs_pinfo
external :: descinit, numroc, pdelset, pdlaprnt, indxl2g, pdgemm
external :: pzelset, pzlaprnt, pdgeadd, pzgesv, pzdotc
integer :: numroc, indxl2g

! Matrix
integer :: Mloc, Nloc, NAkloc, ecc

real(double), allocatable :: Ak(:,:), S(:,:), Yk(:,:)

integer, allocatable :: ipvt(:)
integer :: descAk(9), descMat(9)

integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp, kiniz

real(double) :: t1, t2, dt, alpha, beta

! Initialize a default BLACS context and the processes grid

alpha = 1.0d0
beta = 0.0d0

CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)

Mloc = numroc( total, total/nprow, myrow, 0, nprow )

CALL MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixAk.out',MPI_MODE_RDONLY&
                  &,MPI_INFO_NULL,fa,ierr)

Nloc = numroc( total, total/npcol, mycol, 0, npcol )

call descinit(descMat, total, total, total/nprow, total/nprow,&
             & 0, 0, context, max(1,Mloc), info)

allocate( ipvt( total+(total/nprow) ) )


allocate( S ( Mloc, Nloc ) )

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS.out',MPI_MODE_RDONLY&
                  &,MPI_INFO_NULL,fs,ierr)

DO i = 1, Nloc
   col = INDXL2G(i, descMat(5), mycol, 0, npcol)
!   col = INDXL2G(i, total/npcol, mycol, 0, npcol)
   temp = (col - 1) * total
   IF (descMat(9)>descMat(5)) then
      rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
      offset = (temp + rig -1) * db
      CALL MPI_FILE_READ_AT(fs, offset, S(1,i), descMat(5),&
                           & MPI_DOUBLE_PRECISION,status, ierr)

      rig = INDXL2G(descMat(5)+1, descMat(5), myrow, 0, nprow)
      offset = (temp + rig -1) * db

!      offset = (temp + rig -1 +(descMat(5)*nprow)) * db

      CALL MPI_FILE_READ_AT(fs,offset,S(1+descMat(5),i),descMat(9)-descMat(5),&
                           & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
      offset = (temp+rig-1)*db
      CALL MPI_FILE_READ_AT(fs,offset,S(1,i),descMat(5),&
                           & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO

call MPI_FILE_CLOSE(fs,ierr)

call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixYk.out',MPI_MODE_CREATE+MPI_MODE_WRONLY&
                  &,MPI_INFO_NULL,fy,ierr)

kiniz = 0

DO i = 1, point

   ecc = kdim(i+1) - kdim(i)

   if (ecc <= 0) CYCLE

   ! Computation of local matrix size
   if (ecc < npcol) then
      if (mycol< ecc) then
         NAkloc = 1
      else
         NAkloc = 0
      endif
   else
      NAkloc = numroc( ecc, ecc/npcol, mycol, 0, npcol )
   endif

   allocate( Ak( Mloc, NAkloc ) )
   allocate( Yk( Mloc, NAkloc ) )

   Yk = 0.0d0

   ! Descriptos
   if (ecc < npcol) then
      call descinit(descAk, total, ecc, total/nprow, ecc,&
                   & 0, 0, context, max(1,Mloc), info)
   else
      call descinit(descAk, total, ecc, total/nprow, ecc/npcol,&
                   & 0, 0, context, max(1,Mloc), info)
   endif

   DO j = 1, NAkloc
      if (ecc < npcol) then
         col = INDXL2G(j, ecc, mycol, 0, npcol)
      else
         col = INDXL2G(j, ecc/npcol, mycol, 0, npcol)
      endif
      temp = (col - 1) * total + kiniz
      IF (descAk(9)>descAk(5)) then
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp + rig -1) * db
         CALL MPI_FILE_READ_AT(fa, offset, Ak(1,j), descAk(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
         offset = (temp + rig -1 +(descAk(5)*nprow)) * db
         CALL MPI_FILE_READ_AT(fa,offset,Ak(1+descAk(5),j),descAk(9)-descAk(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
      ELSE
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp+rig-1)*db
         CALL MPI_FILE_READ_AT(fa,offset,Ak(1,j),descAk(5),&
                              & MPI_DOUBLE_PRECISION,status,ierr)
      ENDIF
   ENDDO

   call pdgesv(total, ecc, S, 1, 1, descMat, ipvt, Ak, 1, 1, descAk, info )

   Yk = Ak

write(*,*) 'Yk', Yk

   DO j = 1, NAkloc
      if (ecc < npcol) then
         col = INDXL2G(j, ecc, mycol, 0, npcol)
      else
         col = INDXL2G(j, ecc/npcol, mycol, 0, npcol)
      endif
      temp = (col - 1) * total + kiniz
      IF (descAk(9)>descAk(5)) then
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp + rig -1) * db
         CALL MPI_FILE_WRITE_AT(fy, offset, Yk(1,j), descAk(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
         offset = (temp + rig -1 +(descAk(5)*nprow)) * db
         CALL MPI_FILE_WRITE_AT(fy,offset,Yk(1+descAk(5),j),descAk(9)-descAk(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
      ELSE
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp+rig-1)*db
         CALL MPI_FILE_WRITE_AT(fy,offset,Yk(1,j),descAk(5),&
                              & MPI_DOUBLE_PRECISION,status,ierr)
      ENDIF
   ENDDO

   kiniz = kiniz + ecc*total

   deallocate ( Ak )
   deallocate ( Yk )

enddo

deallocate ( S )
deallocate(ipvt)

call MPI_FILE_CLOSE(fa,ierr)
call MPI_FILE_CLOSE(fy,ierr)

END SUBROUTINE
