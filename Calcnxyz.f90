SUBROUTINE Calcnxyz(context,nprow,npcol,myrow,mycol,total,descVec,&
                   &Mloc,nx,ny,nz,atomtype,numatom,atominv,atomfit,&
                   &totmaxfit,invfit,maxfit,intfit,nn)

use mpi

implicit none

!input
integer, intent(in) :: context,nprow,npcol,myrow,mycol,total
integer, intent(in) :: atomtype,numatom,totmaxfit
integer, intent(in) :: atominv(numatom),atomfit(numatom)
integer, intent(in) :: invfit(atomtype,totmaxfit),maxfit(atomtype)

integer, intent(inout) :: Mloc, descVec(9)
real*8, intent(inout) :: nx(Mloc),ny(Mloc),nz(Mloc),intfit(Mloc),nn

integer :: i, j, m, db,p

! MPI
integer :: ierr, istat
integer :: status(MPI_STATUS_SIZE)
! BLACS/SCALAPACK
integer :: info
external :: blacs_exit, blacs_gridexit, blacs_gridinfo, blacs_get, &
            &blacs_gridinit, blacs_pinfo
external :: descinit, numroc, pdelset, pdlaprnt, indxl2g, pdgemm
external :: pzelset, pzlaprnt, pdgeadd, pzgesv
integer :: numroc, indxl2g

! Matrix

character*4 :: atoms
integer :: ja, atom

real*8 :: valx,valy,valz,valr
complex*16 :: valc

integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp, kiniz

real*8 :: t1, t2, dt, alpha, beta, work( total )


OPEN (UNIT=107, FILE='dipfit', STATUS='OLD', ACTION='READ', IOSTAT=istat)
do i = 1, numatom
   read (107,*) atoms, ja
   atom = atominv(ja)
   do j = 1, maxfit(atom)
      p = invfit(atom,j) 
      READ(107,*) valx, valy, valz
      if ( p==0 ) CYCLE
      call pdelset( nx, p+atomfit(ja),1, descVec, valx)
      call pdelset( ny, p+atomfit(ja),1, descVec, valy)
      call pdelset( nz, p+atomfit(ja),1, descVec, valz)
   enddo
enddo
CLOSE (107)


intfit = 0.0d0
OPEN (UNIT=108, FILE= 'intfit', STATUS='OLD', ACTION='READ', IOSTAT=istat)
do i = 1, numatom
   read (108,*) atoms, ja
   atom = atominv(ja)
   do j = 1, maxfit(atom)
      p = invfit(atom,j)
      READ(108,*) valr
      if ( p==0 ) CYCLE
      call pdelset( intfit, p+atomfit(ja),1, descVec, valr)
   enddo
enddo
CLOSE (108)

call pddot(total,nn,intfit,1,1,descVec,1,intfit,1,1,descVec,1)
CALL MPI_BCAST(nn, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

END SUBROUTINE
