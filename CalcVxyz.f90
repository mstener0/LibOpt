SUBROUTINE CalcVxyz(fa,fd,context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,dipxyz,& 
                   & descV, Mloc, NVloc, Vx, Vy, Vz, Ndiploc,Mdiploc,descdip,dip)

use mpi

implicit none

!input
integer, intent(in) :: fa,fd,context,nprow,npcol,myrow,mycol, total, point
integer, intent(in) :: kdim(point+1), xyz, dipxyz(3)
integer, intent(in) :: descV(9), Mloc, NVloc, Ndiploc, Mdiploc, descdip(9)

real*8, intent(inout) :: dip(Mdiploc,Ndiploc)
real*8, intent(inout) :: Vx(Mloc,NVloc), Vy(Mloc,NVloc), Vz(Mloc,NVloc)


! MPI
integer :: ierr!, fd
integer :: status(MPI_STATUS_SIZE)
! BLACS/SCALAPACK
integer :: info
external :: blacs_get, blacs_pinfo
external :: descinit, numroc, pdelset, pdlaprnt, indxl2g, pdgemv
integer :: numroc, indxl2g

integer :: i, j, db, ecc, Vkidx, totecc, idx

!Matrix
integer :: descAk(9), Nakloc  
real*8, allocatable :: Ak(:,:)  
real*8 :: alpha, beta, work(kdim(point+1))

integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp, kiniz

alpha = 1.0d0
beta = 0.0d0

ierr = 0

totecc = kdim(point+1)

CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)
!acall MPI_FILE_OPEN(MPI_COMM_WORLD,'dipole.out',MPI_MODE_RDONLY&
!                  &,MPI_INFO_NULL,fd,ierr)


DO j = 1, Ndiploc
   if (3 < npcol) then
      col = INDXL2G(j, 1, mycol, 0, npcol)
   else
      col = INDXL2G(j, 3/npcol, mycol, 0, npcol)
   endif

   temp = (col - 1) * totecc

   IF (descdip(9)>descdip(5)) then
      if (totecc < nprow) then
         rig = INDXL2G(1, 1, myrow, 0, nprow)
      else
         rig = INDXL2G(1, totecc/nprow, myrow, 0, nprow)
      endif

      offset = (temp + rig -1) * db

      CALL MPI_FILE_READ_AT(fd, offset, dip(1,j), descdip(5),&
                           & MPI_DOUBLE_PRECISION,status, ierr)

      offset = (temp + rig -1 +(descdip(5)*nprow)) * db

      CALL MPI_FILE_READ_AT(fd,offset,dip(1+descdip(5),j),descdip(9)-descdip(5),&
                           & MPI_DOUBLE_PRECISION,status, ierr)
   ELSE
      if (totecc < nprow) then
         rig = INDXL2G(1, 1, myrow, 0, nprow)
      else
         rig = INDXL2G(1, totecc/nprow, myrow, 0, nprow)
      endif

      offset = (temp+rig-1)*db

      CALL MPI_FILE_READ_AT(fd,offset,dip(1,j),descdip(5),&
                           & MPI_DOUBLE_PRECISION,status,ierr)
   ENDIF
ENDDO

!call MPI_FILE_CLOSE(fd,ierr)

kiniz = 0

DO i = 1, point

   ecc = kdim(i+1) - kdim(i)

   if (ecc <= 0) CYCLE

   idx = i

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

   ! Descriptos
   if (ecc < npcol) then
      call descinit(descAk, total, ecc, total/nprow, 1,&
                   & 0, 0, context, max(1,Mloc), info)
   else
      call descinit(descAk, total, ecc, total/nprow, ecc/npcol,&
                   & 0, 0, context, max(1,Mloc), info)
   endif

   DO j = 1, NAkloc
      if (ecc < npcol) then
         col = INDXL2G(j, 1, mycol, 0, npcol)
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

   kiniz = kiniz + ecc*total

   CALL pdgemv('N',total,ecc,alpha,Ak,1,1,descAk,dip,kdim(i)+1,1,descdip,&
              & 1,beta,Vx,1,i,descV,1)   
   CALL pdgemv('N',total,ecc,alpha,Ak,1,1,descAk,dip,kdim(i)+1,2,descdip,&
              & 1,beta,Vy,1,i,descV,1)
   CALL pdgemv('N',total,ecc,alpha,Ak,1,1,descAk,dip,kdim(i)+1,3,descdip,&
              & 1,beta,Vz,1,i,descV,1)

   deallocate ( Ak )

enddo


END SUBROUTINE

