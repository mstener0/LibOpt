SUBROUTINE StoreA(total,tempi,N,myrank,nprocs,i,fa,block,idxrig,matrixAk) 

implicit none

include 'mpif.h'

!input

integer, intent(in) :: myrank, nprocs, total, tempi, N, i, fa
integer, intent(in) :: block(N+1), idxrig(0:nprocs)
real*8, intent(in) :: MatrixAk(tempi,block(i+1)-block(i))

! MPI
integer :: ierr, istat

integer :: j, p, l, m, dimcol, tempcol, db

integer, allocatable :: idxcol(:)

real*8, allocatable :: Pr(:,:)


integer :: status(MPI_STATUS_SIZE)
integer (KIND=MPI_OFFSET_KIND) :: offset, temp, tot

real*8 :: t1, t2, dt

CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)

tot = total

if (block(i+1) - block(i) > nprocs) then

   dimcol = block(i+1) - block(i)
   tempcol = dimcol / nprocs
   allocate (idxcol(0:nprocs))
   idxcol = 0
   do l = 1,nprocs-1
      idxcol(l) = idxcol(l-1) + tempcol
   enddo
   idxcol(nprocs) = dimcol
   tempcol = idxcol(myrank+1) - idxcol(myrank)
   allocate (Pr(total, tempcol))
   Pr = 0.

   do l = 0, nprocs-1
      if (myrank==l) then
         do j = 0, nprocs-1
            if (j==l) then
               do p = 1, tempcol
                  do m = 1, tempi
                     Pr(idxrig(j)+m,p) = MatrixAk(m,idxcol(j)+p)
                  enddo
               enddo
            else
               do p = idxcol(j)+1, idxcol(j+1)
                  call MPI_SEND(MatrixAk(1,p),&
                  & tempi,MPI_DOUBLE_PRECISION,j, &
                  & myrank,MPI_COMM_WORLD,ierr)
               enddo
            endif
         enddo
      else
         do p = 1, tempcol
            call MPI_RECV(Pr(idxrig(l)+1,p),&
            & idxrig(l+1)-idxrig(l),MPI_DOUBLE_PRECISION,l, &
            & l,MPI_COMM_WORLD,status,ierr)
         enddo
      endif
   enddo

   offset = (block(i)+idxcol(myrank))*tot*db
   CALL MPI_FILE_WRITE_AT(fa, offset, Pr(1,1), total*tempcol,&
   & MPI_DOUBLE_PRECISION,status, ierr)

   deallocate (Pr)
   deallocate (idxcol)

else

   DO j = block(i)+1, block(i+1)
      temp = (j - 1) * tot
      offset = (temp + idxrig(myrank)) * db
      CALL MPI_FILE_WRITE_AT_ALL(fa, offset, MatrixAk(1,j-block(i)), tempi,&
                                & MPI_DOUBLE_PRECISION,status, ierr)
   ENDDO

endif

END SUBROUTINE StoreA
