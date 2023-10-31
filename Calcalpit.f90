SUBROUTINE Calcalpit(context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,&
          &dipxyz,atomtype,numatom,atominv,atomfit,totmaxfit,   &
          &invfit,maxfit,intervals,Ein,Efin,numstep,wi,YN,numpoint, &
          &enpoint)


use mpi

implicit none

!input
integer, intent(in) :: context,nprow,npcol,myrow,mycol, total, point
integer, intent(in) :: kdim(point+1), xyz, dipxyz(3)
integer, intent(in) :: atomtype,numatom,atominv(numatom)
integer, intent(in) :: atomfit(numatom),totmaxfit, maxfit(atomtype)
integer, intent(in) :: invfit(atomtype,totmaxfit)
real*8, intent(in) :: intervals(point+1)
real*8, intent(inout) :: Ein, Efin, wi, enpoint(numpoint)
integer, intent(inout) :: numstep
integer, intent(in) :: numpoint
character*1, intent(in) :: YN

integer, parameter :: double = selected_real_kind(14)
real(double), parameter :: eVha = 27.2113961
real(double), parameter :: epsi = 1.0d-6

integer :: i, j, m, db, ien, k, istat, ite

! MPI
integer :: ierr, fa, fh, fs, fy
integer :: status(MPI_STATUS_SIZE)
! BLACS/SCALAPACK
integer :: info
external :: blacs_exit, blacs_gridexit, blacs_gridinfo, blacs_get, &
            &blacs_gridinit, blacs_pinfo
external :: descinit, numroc, pdelset, pdlaprnt, indxl2g, pdgemm
external :: pzelset, pzlaprnt, pdgeadd, pzgesv, pzdotc
integer :: numroc, indxl2g

! Matrix
integer :: Mloc, Nloc, NAkloc, ecc, NVloc, Ndiploc, Mdiploc

real(double), allocatable :: Ak(:,:), Dk(:,:), Qr(:,:), Qi(:,:), L(:,:), work(:)
real(double), allocatable :: Mr(:,:), Mi(:,:), S(:,:), Vx(:,:), Vy(:,:), Vz(:,:)
real(double), allocatable :: dxr(:), dxi(:), dyr(:), dyi(:), dzr(:), dzi(:)
real(double), allocatable :: nx(:), ny(:), nz(:), wr(:)
real(double), allocatable :: axr(:), axi(:), ayr(:), ayi(:), azr(:), azi(:)
real(double), allocatable :: bxr(:), bxi(:), byr(:), byi(:), bzr(:), bzi(:)
real(double), allocatable :: vectr(:), vecti(:), dip(:,:), vecAr(:), vecAi(:)
real(double), allocatable :: intfit(:), testr(:), testi(:)
real(double), allocatable :: Yk(:,:)
real(double), allocatable :: bxrp(:), bxip(:), byrp(:), byip(:), bzrp(:), bzip(:)

integer, allocatable :: ipvt(:)
integer :: descAk(9), descMat(9), descV(9)
integer :: descVec(9), descdip(9), descvect(9)
integer :: descVinc(9)

real(double) :: aval, val, skr, ski, valr, vali, Ek
real(double) :: step, nn, deltar, deltai, deltartot, deltaitot
integer :: totecc

complex(double) :: sk, w, term1, term2, valc, valc2, tempc1, tempc2
complex(double), allocatable :: C(:,:), dx(:), dy(:), dz(:)
complex(double), allocatable :: ax(:), ay(:), az(:)


complex(double), allocatable :: bx(:), by(:), bz(:)
complex(double) :: valalp,valbet
complex(double), allocatable :: intfitcompl(:)
complex(double), allocatable :: dxnxvinc(:,:),dynyvinc(:,:),dznzvinc(:,:)

integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp, kiniz

real(double) :: t1, t2, dt, alpha, beta, t3, t4, dt2

! Initialize a default BLACS context and the processes grid

t3 = MPI_wtime()

alpha = 1.0d0
beta = 0.0d0

!OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!& POSITION='APPEND',IOSTAT=istat)

CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, db, ierr)

Mloc = numroc( total, total/nprow, myrow, 0, nprow )

if (point < npcol) then
   if (mycol< point) then
      NVloc = 1
   else
      NVloc = 0
   endif
else
   NVloc = numroc( point, point/npcol, mycol, 0, npcol )
endif

allocate( Vx ( Mloc, NVloc ) )
allocate( Vy ( Mloc, NVloc ) )
allocate( Vz ( Mloc, NVloc ) )

if (point < npcol) then
   call descinit(descV, total, point, total/nprow, point,&
             & 0, 0, context, max(1,Mloc), info)
else
   call descinit(descV, total, point, total/nprow, point/npcol,&
                & 0, 0, context, max(1,Mloc), info)
endif

Vx = 0.0d0
Vy = 0.0d0
Vz = 0.0d0

CALL MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixAk.out',MPI_MODE_RDONLY&
                  &,MPI_INFO_NULL,fa,ierr)

!write(*,*) 'kdim', kdim

if (3 < npcol) then
   if (mycol < 3) then
      Ndiploc = 1
   else
      Ndiploc = 0
   endif
else
   Ndiploc = numroc( 3, 3/npcol, mycol, 0, npcol )
endif


totecc = kdim(point+1)
if (totecc < nprow) then
   if (myrow < totecc) then
      Mdiploc = 1
   else
      Mdiploc = 0
   endif
else
   Mdiploc = numroc( totecc, totecc/nprow, myrow, 0, nprow )
endif


allocate ( dip ( Mdiploc, Ndiploc ) )
allocate ( vecAr ( Mdiploc ) )
allocate ( vecAi ( Mdiploc ) )
allocate ( testr ( Mdiploc ) )
allocate ( testi ( Mdiploc ) )



if ( totecc < nprow) then
   call descinit(descvect, totecc, 1, totecc, 1,&
                & 0, 0, context, max(1,Mdiploc), info)
else
   call descinit(descvect, totecc, 1, totecc/nprow, 1,&
                & 0, 0, context, max(1,Mdiploc), info)
endif




if (3 < npcol) then
   if ( totecc < nprow) then
      call descinit(descdip, totecc, 3, totecc, 3,&
                   & 0, 0, context, max(1,Mdiploc), info)
   else
      call descinit(descdip, totecc, 3, totecc/nprow, 3,&
                   & 0, 0, context, max(1,Mdiploc), info)
   endif
else
   if ( totecc < nprow) then
      call descinit(descdip, totecc, 3, totecc, 3/npcol,&
                   & 0, 0, context, max(1,Mdiploc), info)
   else
      call descinit(descdip, totecc, 3, totecc/nprow, 3/npcol,&
                   & 0, 0, context, max(1,Mdiploc), info)
   endif
endif



CALL CalcVxyz(fa,context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,dipxyz,&
          & descV, Mloc, NVloc, Vx, Vy, Vz, Ndiploc, Mdiploc, descdip, dip)

Nloc = numroc( total, total/npcol, mycol, 0, npcol )

call descinit(descMat, total, total, total/nprow, total/nprow,&
             & 0, 0, context, max(1,Mloc), info)
call descinit(descVec, total, 1, total/nprow, 1, 0, 0, context&
             &, max(1,Mloc), info)
call descinit(descVinc, total, 2, total/nprow, 2, 0, 0, context&
             &, max(1,Mloc), info)

if (YN == 'Y') then
allocate ( intfitcompl ( Mloc ) )
intfitcompl = 0.0d0
endif


allocate ( nx ( Mloc ) )
allocate ( ny ( Mloc ) )
allocate ( nz ( Mloc ) )

allocate ( intfit ( Mloc ) ) 

nx = 0.0d0
ny = 0.0d0
nz = 0.0d0


CALL Calcnxyz(context,nprow,npcol,myrow,mycol,total,descVec,Mloc,nx,ny,nz,&
             &atomtype,numatom,atominv,atomfit,totmaxfit,invfit,maxfit,&
             &intfit,nn)

if (YN == 'Y') then
intfitcompl = CMPLX(intfit,0.0d0)
endif

Call CalcYk(context,nprow,npcol,myrow,mycol,total,point,kdim)

CALL MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixYk.out',MPI_MODE_RDONLY&
                  &,MPI_INFO_NULL,fy,ierr)


!inserire ciclo energie

Ein = Ein/eVha
Efin = Efin/Evha
allocate (wr(numstep+1))
wr(1) = Ein
wr(numstep+1) = Efin
step = (Efin-Ein)/numstep
do i = 2, numstep
   wr(i) = wr(i-1) + step
enddo

wi = wi/Evha

allocate ( ax ( numstep ) )
allocate ( ay ( numstep ) )
allocate ( az ( numstep ) )
allocate ( axr ( numstep ) )
allocate ( axi ( numstep ) )
allocate ( ayr ( numstep ) )
allocate ( ayi ( numstep ) )
allocate ( azr ( numstep ) )
allocate ( azi ( numstep ) )
allocate( ipvt( total+(total/nprow) ) )

axr = 0.0d0
axi = 0.0d0
ax = 0.0d0
ayr = 0.0d0
ayi = 0.0d0
ay = 0.0d0
azr = 0.0d0
azi = 0.0d0
az = 0.0d0

DO ien = 1, numstep
      
!   allocate( Qr ( Mloc, Nloc ) )
!   allocate( Qi ( Mloc, Nloc ) )
!   allocate( Dk ( Mloc, Nloc ) )
   allocate( work( total ) )
   allocate( dxr ( Mloc ) )
   allocate( dxi ( Mloc ) )
   allocate( dyr ( Mloc ) )
   allocate( dyi ( Mloc ) )
   allocate( dzr ( Mloc ) )
   allocate( dzi ( Mloc ) )


!   Qr = 0.0d0
!   Qi = 0.0d0
   dxr = 0.0d0
   dxi = 0.0d0
   dyr = 0.0d0
   dyi = 0.0d0
   dzr = 0.0d0
   dzi = 0.0d0

   w = CMPLX(wr(ien),wi)

!   kiniz = 0

   DO i = 1, point
   
      ecc = kdim(i+1) - kdim(i)
   
      if (ecc <= 0) CYCLE
      Ek = (intervals(i) + intervals(i+1))/2
      tempc1 = 1.0d0
      tempc2 = Ek
      term1 = tempc1/(tempc2-w)
      term2 = tempc1/(tempc2+w)
      sk = 2*(term1+term2)

      skr = REAL(sk)
      ski = AIMAG(sk)


!      ! Computation of local matrix size
!      if (ecc < npcol) then
!         if (mycol< ecc) then
!            NAkloc = 1
!         else
!            NAkloc = 0
!         endif
!      else
!         NAkloc = numroc( ecc, ecc/npcol, mycol, 0, npcol )
!      endif 

!      allocate( Ak( Mloc, NAkloc ) )
!
!      ! Descriptos
!      if (ecc < npcol) then
!         call descinit(descAk, total, ecc, total/nprow, ecc,&
!                      & 0, 0, context, max(1,Mloc), info)
!      else
!         call descinit(descAk, total, ecc, total/nprow, ecc/npcol,&
!                      & 0, 0, context, max(1,Mloc), info)
!      endif
!
!      DO j = 1, NAkloc
!         if (ecc < npcol) then  
!            col = INDXL2G(j, ecc, mycol, 0, npcol)
!         else      
!            col = INDXL2G(j, ecc/npcol, mycol, 0, npcol)
!         endif
!         temp = (col - 1) * total + kiniz
!         IF (descAk(9)>descAk(5)) then
!            rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
!            offset = (temp + rig -1) * db
!            CALL MPI_FILE_READ_AT(fa, offset, Ak(1,j), descAk(5),&
!                                 & MPI_DOUBLE_PRECISION,status, ierr)
!            offset = (temp + rig -1 +(descAk(5)*nprow)) * db
!            CALL MPI_FILE_READ_AT(fa,offset,Ak(1+descAk(5),j),descAk(9)-descAk(5),&
!                                 & MPI_DOUBLE_PRECISION,status, ierr)
!         ELSE
!            rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
!            offset = (temp+rig-1)*db
!            CALL MPI_FILE_READ_AT(fa,offset,Ak(1,j),descAk(5),&
!                                 & MPI_DOUBLE_PRECISION,status,ierr)
!         ENDIF
!      ENDDO
!
!      kiniz = kiniz + ecc*total

!      CALL pdgemm('N','T',total,total,ecc,alpha,Ak,1,1,descAk,Ak,1,1,descAk,&
!                 & beta,Dk,1,1,descMat)

!      deallocate ( Ak )

!      call pdgeadd('N', total,total,skr,Dk,1,1,descMat,alpha,Qr,1,1,descMat)
!      call pdgeadd('N', total,total,ski,Dk,1,1,descMat,alpha,Qi,1,1,descMat)

      call pdgeadd('N', total,1,skr,Vx,1,i,descV,alpha,dxr,1,1,descVec)
      call pdgeadd('N', total,1,ski,Vx,1,i,descV,alpha,dxi,1,1,descVec)
      call pdgeadd('N', total,1,skr,Vy,1,i,descV,alpha,dyr,1,1,descVec)
      call pdgeadd('N', total,1,ski,Vy,1,i,descV,alpha,dyi,1,1,descVec)
      call pdgeadd('N', total,1,skr,Vz,1,i,descV,alpha,dzr,1,1,descVec)
      call pdgeadd('N', total,1,ski,Vz,1,i,descV,alpha,dzi,1,1,descVec)

   enddo

   allocate ( dx ( Mloc ) )
   allocate ( dy ( Mloc ) )
   allocate ( dz ( Mloc ) )

   dx = CMPLX (dxr, dxi )
   dy = CMPLX (dyr, dyi )
   dz = CMPLX (dzr, dzi )


   if (YN == 'Y') THEN
      allocate ( dxnxvinc ( Mloc,2 ) )
      allocate ( dynyvinc ( Mloc,2 ) )
      allocate ( dznzvinc ( Mloc,2 ) )

      dxnxvinc(:,1) = dx
      dynyvinc(:,1) = dy
      dznzvinc(:,1) = dz

      dxnxvinc(:,2) = CMPLX(intfit, 0.0d0)
      dynyvinc(:,2) = CMPLX(intfit, 0.0d0)
      dznzvinc(:,2) = CMPLX(intfit, 0.0d0)
   endif

   deallocate(dxr,dxi,dyr,dyi,dzr,dzi)
!   deallocate(Dk)

   allocate( L ( Mloc, Nloc ) )

   call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixL.out',MPI_MODE_RDONLY&
                     &,MPI_INFO_NULL,fh,ierr)

   DO i = 1, Nloc
      col = INDXL2G(i, total/nprow, mycol, 0, npcol)
      temp = (col - 1) * total
      IF (descMat(9)>descMat(5)) then
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp + rig -1) * db
         CALL MPI_FILE_READ_AT(fh, offset, L(1,i), descMat(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
         offset = (temp + rig -1 +(descMat(5)*nprow)) * db
         CALL MPI_FILE_READ_AT(fh,offset,L(1+descMat(5),i),descMat(9)-descMat(5),&
                               & MPI_DOUBLE_PRECISION,status, ierr)
      ELSE
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp+rig-1)*db
         CALL MPI_FILE_READ_AT(fh,offset,L(1,i),descMat(5),&
                              & MPI_DOUBLE_PRECISION,status,ierr)
      ENDIF
   ENDDO

   call MPI_FILE_CLOSE(fh,ierr)

!   allocate( Mr ( Mloc, Nloc ) )
!   allocate( Mi ( Mloc, Nloc ) )
!
!   CALL pdgemm('N','N',total,total,total,alpha,Qr,1,1,descMat,L,1,1,descMat,&
!              & beta,Mr,1,1,descMat)
!
!   CALL pdgemm('N','N',total,total,total,alpha,Qi,1,1,descMat,L,1,1,descMat,&
!              & beta,Mi,1,1,descMat)
!
!
!   deallocate(Qr, Qi)
!
!   allocate( S ( Mloc, Nloc ) )
!
!   call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS.out',MPI_MODE_RDONLY&
!                     &,MPI_INFO_NULL,fs,ierr)
!
!   DO i = 1, Nloc
!      col = INDXL2G(i, total/nprow, mycol, 0, npcol)
!      temp = (col - 1) * total
!      IF (descMat(9)>descMat(5)) then
!         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
!         offset = (temp + rig -1) * db
!         CALL MPI_FILE_READ_AT(fs, offset, S(1,i), descMat(5),&
!                              & MPI_DOUBLE_PRECISION,status, ierr)
!         offset = (temp + rig -1 +(descMat(5)*nprow)) * db
!         CALL MPI_FILE_READ_AT(fs,offset,S(1+descMat(5),i),descMat(9)-descMat(5),&
!                              & MPI_DOUBLE_PRECISION,status, ierr)
!      ELSE
!         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
!         offset = (temp+rig-1)*db
!         CALL MPI_FILE_READ_AT(fs,offset,S(1,i),descMat(5),&
!                              & MPI_DOUBLE_PRECISION,status,ierr)
!      ENDIF
!   ENDDO
!
!   call MPI_FILE_CLOSE(fs,ierr)
!
!   allocate( C ( Mloc, Nloc ) )
!
!   do j = 1, Nloc
!      do m = 1, Mloc
!         valr = S(m,j)   +Mr(m,j)        !  -Mr(m,j)
!         vali = 0.0d0    +Mi(m,j)        !  -Mi(m,j)
!         C(m,j) = CMPLX ( valr , vali )
!      enddo
!   enddo
!
!   deallocate(Mr,Mi,S)

   t1 = MPI_wtime()

   do m = 1, xyz
      if (dipxyz(m) == 1) then


         allocate (bxr ( Mloc ) )
         allocate (bxi ( Mloc ) )

         allocate (bxrp ( Mloc ) )
         allocate (bxip ( Mloc ) ) 

         bxr = 0.0d0
         bxi = 0.0d0

         bxrp = 0.0d0
         bxip = 0.0d0

! CICLO ITERAZIONI
            
         allocate (vectr ( Mloc ) )
         allocate (vecti ( Mloc ) )

         do ite = 1, 100

            CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bxr,1,1,descVec,&
                       & 1,beta,vectr,1,1,descVec,1)
            CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bxi,1,1,descVec,&
                       & 1,beta,vecti,1,1,descVec,1)


            kiniz = 0

            DO i = 1, point

               ecc = kdim(i+1) - kdim(i)
               if (ecc <= 0) CYCLE

               Ek = (intervals(i) + intervals(i+1))/2
               tempc1 = CMPLX(1.0d0, 0.0d0)
               tempc2 = Ek - w
               term1 = tempc1/tempc2

               tempc2 = Ek
               term2 = tempc1/(tempc2+w)
               sk = 2*(term1+term2)

               skr = REAL(sk)
               ski = AIMAG(sk)


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
                     CALL MPI_FILE_READ_AT(fa,offset,Ak(1+descAk(5),j),&
                                          &descAk(9)-descAk(5),MPI_DOUBLE_PRECISION,&
                                          &status, ierr)
                  ELSE
                     rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                     offset = (temp+rig-1)*db
                     CALL MPI_FILE_READ_AT(fa,offset,Ak(1,j),descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status,ierr)
                  ENDIF
               ENDDO

               CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectr,1,1,descVec,&
                          & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
               CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vecti,1,1,descVec,&
                          & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

               vecAr = - vecAr
               vecAi = - vecAi

               call pdgeadd('N',ecc,1,alpha,dip,kdim(i)+1,1,descdip,alpha,vecAr,kdim(i)+1,1,descvect)

               testr = skr * vecAr
               testr = testr - ( ski* vecAi)

               testi = skr * vecAi
               testi = testi + (ski * vecAr)
 
               vecAr = testr
               vecAi = testi

               deallocate ( Ak )

               allocate( Yk( Mloc, NAkloc ) )

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
                     CALL MPI_FILE_READ_AT(fy, offset, Yk(1,j), descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status, ierr)
                     offset = (temp + rig -1 +(descAk(5)*nprow)) * db
                     CALL MPI_FILE_READ_AT(fy,offset,Yk(1+descAk(5),j),&
                                          &descAk(9)-descAk(5),MPI_DOUBLE_PRECISION,&
                                          &status, ierr)
                  ELSE
                     rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                     offset = (temp+rig-1)*db
                     CALL MPI_FILE_READ_AT(fy,offset,Yk(1,j),descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status,ierr)
                  ENDIF
               ENDDO

               kiniz = kiniz + ecc*total

               CALL pdgemv('N',total,ecc,alpha,Yk,1,1,descAk,vecAr,kdim(i)+1,1,descvect,&
                          & 1,beta,bxrp,1,1,descVec,1)
               CALL pdgemv('N',total,ecc,alpha,Yk,1,1,descAk,vecAi,kdim(i)+1,1,descvect,&
                          & 1,beta,bxip,1,1,descVec,1)

               deallocate ( Yk )

            enddo

            deltar = 0.0d0
            deltai = 0.0d0

            do j = 1, Mloc

               deltar = deltar + ABS(bxrp(j) - bxr(j))
               deltai = deltai + ABS(bxip(j) - bxi(j))

            enddo

CALL MPI_ALLREDUCE(deltar,deltartot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(deltai,deltaitot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

            bxr = 0.2*bxrp + 0.8*bxr
            bxi = 0.2*bxip + 0.8*bxi

            IF ( deltartot < epsi .and. deltaitot < epsi) exit

            bxrp = 0.0d0
            bxip = 0.0d0

         enddo

         deallocate(bxrp,bxip)

         call pddot(totecc,valr,dip,1,1,descdip,1,vecAr,1,1,descvect,1)
         call pddot(totecc,vali,dip,1,1,descdip,1,vecAi,1,1,descvect,1)

         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(vali, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)

         axr(ien) = axr(ien) + valr      
         axi(ien) = axi(ien) + vali     

         do i = 1, numpoint
            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
               do j = 0, nprow-1
                  if (myrow == j .and. mycol ==0) then
OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)
                     if (j == 0) write(200,*) 'point', i, 'energy', enpoint(i)
                     if (j == 0) write(200,*) 'bx'
                     do k = 1, Mloc
                        write(200,*) bxr(k), bxi(k)
                     enddo
CLOSE(200)
                  endif
                  call MPI_Barrier( mpi_comm_world, ierr)
               enddo
            endif
         enddo

         do i = 1, numpoint
            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
               do j = 0, nprow-1 
                  if (myrow == j .and. mycol ==0) then
OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)
                     if (j == 0) write(200,*) 'Rirx'
                     do k = 1, Mdiploc
                        write(200,*) vecAr(k), vecAi(k)
                     enddo
CLOSE(200)
                  endif
                  call MPI_Barrier( mpi_comm_world, ierr)
               enddo
            endif
         enddo
 
         deallocate (vectr,vecti)

         call pddot(total,valr,bxr,1,1,descVec,1,nx,1,1,descVec,1)
         call pddot(total,vali,bxi,1,1,descVec,1,nx,1,1,descVec,1)

         deallocate(bxr,bxi)

         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         ax(ien) = ax(ien) + CMPLX(valr, vali)


      else if (dipxyz(m) == 2 ) then
 
         allocate (byr ( Mloc ) )
         allocate (byi ( Mloc ) )

         allocate (byrp ( Mloc ) )
         allocate (byip ( Mloc ) )

         byr = 0.0d0
         byi = 0.0d0

         byrp = 0.0d0
         byip = 0.0d0

! CICLO ITERAZIONI

         allocate (vectr ( Mloc ) )
         allocate (vecti ( Mloc ) )

         do ite = 1, 100

            CALL pdgemv('N',total,total,alpha,L,1,1,descMat,byr,1,1,descVec,&
                       & 1,beta,vectr,1,1,descVec,1)
            CALL pdgemv('N',total,total,alpha,L,1,1,descMat,byi,1,1,descVec,&
                       & 1,beta,vecti,1,1,descVec,1)


            kiniz = 0

            DO i = 1, point

               ecc = kdim(i+1) - kdim(i)
               if (ecc <= 0) CYCLE

               Ek = (intervals(i) + intervals(i+1))/2
               tempc1 = CMPLX(1.0d0, 0.0d0)
               tempc2 = Ek - w
               term1 = tempc1/tempc2

               tempc2 = Ek
               term2 = tempc1/(tempc2+w)
               sk = 2*(term1+term2)

               skr = REAL(sk)
               ski = AIMAG(sk)


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
                     CALL MPI_FILE_READ_AT(fa,offset,Ak(1+descAk(5),j),&
                                          &descAk(9)-descAk(5),MPI_DOUBLE_PRECISION,&
                                          &status, ierr)
                  ELSE
                     rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                     offset = (temp+rig-1)*db
                     CALL MPI_FILE_READ_AT(fa,offset,Ak(1,j),descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status,ierr)
                  ENDIF
               ENDDO

               CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectr,1,1,descVec,&
                          & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
               CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vecti,1,1,descVec,&
                          & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

               vecAr = - vecAr
               vecAi = - vecAi

               call pdgeadd('N',ecc,1,alpha,dip,kdim(i)+1,2,descdip,alpha,vecAr,kdim(i)+1,1,descvect)

               testr = skr * vecAr
               testr = testr - ( ski* vecAi)

               testi = skr * vecAi
               testi = testi + (ski * vecAr)

               vecAr = testr
               vecAi = testi

               deallocate ( Ak )

               allocate( Yk( Mloc, NAkloc ) )

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
                     CALL MPI_FILE_READ_AT(fy, offset, Yk(1,j), descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status, ierr)
                     offset = (temp + rig -1 +(descAk(5)*nprow)) * db
                     CALL MPI_FILE_READ_AT(fy,offset,Yk(1+descAk(5),j),&
                                          &descAk(9)-descAk(5),MPI_DOUBLE_PRECISION,&
                                          &status, ierr)
                  ELSE
                     rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                     offset = (temp+rig-1)*db
                     CALL MPI_FILE_READ_AT(fy,offset,Yk(1,j),descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status,ierr)
                  ENDIF
               ENDDO

               kiniz = kiniz + ecc*total

               CALL pdgemv('N',total,ecc,alpha,Yk,1,1,descAk,vecAr,kdim(i)+1,1,descvect,&
                          & 1,beta,byrp,1,1,descVec,1)
               CALL pdgemv('N',total,ecc,alpha,Yk,1,1,descAk,vecAi,kdim(i)+1,1,descvect,&
                          & 1,beta,byip,1,1,descVec,1)

               deallocate ( Yk )

            enddo

            deltar = 0.0d0
            deltai = 0.0d0

            do j = 1, Mloc

               deltar = deltar + ABS(byrp(j) - byr(j))
               deltai = deltai + ABS(byip(j) - byi(j))

            enddo

CALL MPI_ALLREDUCE(deltar,deltartot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(deltai,deltaitot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

            byr = 0.2*byrp + 0.8*byr
            byi = 0.2*byip + 0.8*byi

            IF ( deltartot < epsi .and. deltaitot < epsi) exit

            byrp = 0.0d0
            byip = 0.0d0

         enddo

         deallocate (byrp, byip)

         call pddot(totecc,valr,dip,1,2,descdip,1,vecAr,1,1,descvect,1)
         call pddot(totecc,vali,dip,1,2,descdip,1,vecAi,1,1,descvect,1)

         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(vali, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)

         ayr(ien) = ayr(ien) + valr
         ayi(ien) = ayi(ien) + vali

         do i = 1, numpoint
            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
               do j = 0, nprow-1
                  if (myrow == j .and. mycol ==0) then
OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)
                     if (j == 0) write(200,*) 'point', i, 'energy', enpoint(i)
                     if (j == 0) write(200,*) 'by'
                     do k = 1, Mloc
                        write(200,*) byr(k), byi(k)
                     enddo
CLOSE(200)
                  endif
                  call MPI_Barrier( mpi_comm_world, ierr)
               enddo
            endif
         enddo

         do i = 1, numpoint
            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
               do j = 0, nprow-1
                  if (myrow == j .and. mycol ==0) then
OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)
                     if (j == 0) write(200,*) 'Riry'
                     do k = 1, Mdiploc
                        write(200,*) vecAr(k), vecAi(k)
                     enddo
CLOSE(200)
                  endif
                  call MPI_Barrier( mpi_comm_world, ierr)
               enddo
            endif
         enddo

         deallocate (vectr,vecti)

         call pddot(total,valr,byr,1,1,descVec,1,ny,1,1,descVec,1)
         call pddot(total,vali,byi,1,1,descVec,1,ny,1,1,descVec,1)

         deallocate(byr,byi)

         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         ay(ien) = ay(ien) + CMPLX(valr, vali)



      else if (dipxyz(m) == 3 ) then

         allocate (bzr ( Mloc ) )
         allocate (bzi ( Mloc ) )

         allocate (bzrp ( Mloc ) )
         allocate (bzip ( Mloc ) )

         bzr = 0.0d0
         bzi = 0.0d0

         bzrp = 0.0d0
         bzip = 0.0d0

! CICLO ITERAZIONI

         allocate (vectr ( Mloc ) )
         allocate (vecti ( Mloc ) )

         do ite = 1, 100

            CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bzr,1,1,descVec,&
                       & 1,beta,vectr,1,1,descVec,1)
            CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bzi,1,1,descVec,&
                       & 1,beta,vecti,1,1,descVec,1)


            kiniz = 0

            DO i = 1, point

               ecc = kdim(i+1) - kdim(i)
               if (ecc <= 0) CYCLE

               Ek = (intervals(i) + intervals(i+1))/2
               tempc1 = CMPLX(1.0d0, 0.0d0)
               tempc2 = Ek - w
               term1 = tempc1/tempc2

               tempc2 = Ek
               term2 = tempc1/(tempc2+w)
               sk = 2*(term1+term2)

               skr = REAL(sk)
               ski = AIMAG(sk)


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
                     CALL MPI_FILE_READ_AT(fa,offset,Ak(1+descAk(5),j),&
                                          &descAk(9)-descAk(5),MPI_DOUBLE_PRECISION,&
                                          &status, ierr)
                  ELSE
                     rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                     offset = (temp+rig-1)*db
                     CALL MPI_FILE_READ_AT(fa,offset,Ak(1,j),descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status,ierr)
                  ENDIF
               ENDDO

               CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectr,1,1,descVec,&
                          & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
               CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vecti,1,1,descVec,&
                          & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

               vecAr = - vecAr
               vecAi = - vecAi

               call pdgeadd('N',ecc,1,alpha,dip,kdim(i)+1,3,descdip,alpha,vecAr,kdim(i)+1,1,descvect)

               testr = skr * vecAr
               testr = testr - ( ski* vecAi)

               testi = skr * vecAi
               testi = testi + (ski * vecAr)

               vecAr = testr
               vecAi = testi

               deallocate ( Ak )

               allocate( Yk( Mloc, NAkloc ) )

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
                     CALL MPI_FILE_READ_AT(fy, offset, Yk(1,j), descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status, ierr)
                     offset = (temp + rig -1 +(descAk(5)*nprow)) * db
                     CALL MPI_FILE_READ_AT(fy,offset,Yk(1+descAk(5),j),&
                                          &descAk(9)-descAk(5),MPI_DOUBLE_PRECISION,&
                                          &status, ierr)
                  ELSE
                     rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                     offset = (temp+rig-1)*db
                     CALL MPI_FILE_READ_AT(fy,offset,Yk(1,j),descAk(5),&
                                          & MPI_DOUBLE_PRECISION,status,ierr)
                  ENDIF
               ENDDO

               kiniz = kiniz + ecc*total

               CALL pdgemv('N',total,ecc,alpha,Yk,1,1,descAk,vecAr,kdim(i)+1,1,descvect,&
                          & 1,beta,bzrp,1,1,descVec,1)
               CALL pdgemv('N',total,ecc,alpha,Yk,1,1,descAk,vecAi,kdim(i)+1,1,descvect,&
                          & 1,beta,bzip,1,1,descVec,1)

               deallocate ( Yk )

            enddo

            deltar = 0.0d0
            deltai = 0.0d0

            do j = 1, Mloc

               deltar = deltar + ABS(bzrp(j) - bzr(j))
               deltai = deltai + ABS(bzip(j) - bzi(j))

            enddo

CALL MPI_ALLREDUCE(deltar,deltartot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(deltai,deltaitot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

write(*,*) 'deltar', deltar, 'deltai', deltai


            bzr = 0.1*bzrp + 0.9*bzr
            bzi = 0.1*bzip + 0.9*bzi

            IF ( deltartot < epsi .and. deltaitot < epsi) exit

            bzrp = 0.0d0
            bzip = 0.0d0

         enddo

         deallocate (bzrp, bzip)

         call pddot(totecc,valr,dip,1,3,descdip,1,vecAr,1,1,descvect,1)
         call pddot(totecc,vali,dip,1,3,descdip,1,vecAi,1,1,descvect,1)

         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(vali, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)

         azr(ien) = azr(ien) + valr
         azi(ien) = azi(ien) + vali

         do i = 1, numpoint
            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
               do j = 0, nprow-1
                  if (myrow == j .and. mycol ==0) then
OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)
                     if (j == 0) write(200,*) 'point', i, 'energy', enpoint(i)
                     if (j == 0) write(200,*) 'bz'
                     do k = 1, Mloc
                        write(200,*) bzr(k), bzi(k)
                     enddo
CLOSE(200)
                  endif
                  call MPI_Barrier( mpi_comm_world, ierr)
               enddo
            endif
         enddo

         do i = 1, numpoint
            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
               do j = 0, nprow-1
                  if (myrow == j .and. mycol ==0) then
OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
& POSITION='APPEND',IOSTAT=istat)
                     if (j == 0) write(200,*) 'Rirz'
                     do k = 1, Mdiploc
                        write(200,*) vecAr(k), vecAi(k)
                     enddo
CLOSE(200)
                  endif
                  call MPI_Barrier( mpi_comm_world, ierr)
               enddo
            endif
         enddo

         deallocate (vectr,vecti)

         call pddot(total,valr,bzr,1,1,descVec,1,nz,1,1,descVec,1)
         call pddot(total,vali,bzi,1,1,descVec,1,nz,1,1,descVec,1)

         deallocate(bzr,bzi)

         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(valr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         az(ien) = az(ien) + CMPLX(valr, vali)


      endif
   enddo

   t2 = MPI_wtime()
   dt = t2 -t1

   deallocate(L)
   deallocate(dx,dy,dz)
!   deallocate( C )
   deallocate( work )

if (YN=='Y') then
   deallocate(dxnxvinc,dynyvinc,dznzvinc)
endif

enddo

!write(*,*) 'dip', dip

deallocate(dip,vecAr,vecAi,testr, testi)
deallocate(nx, ny, nz)
if(YN=='Y') then
deallocate(intfitcompl)
endif

deallocate(ipvt)

deallocate(Vx,Vy,Vz)

!chiudere ciclo energie

call MPI_FILE_CLOSE(fa,ierr)
call MPI_FILE_CLOSE(fy,ierr)

if (myrow==0 .and. mycol==0) then
   do m = 1, xyz
      if (dipxyz(m) == 1) then
         write(*,*) 'ax'
         do i = 1, numstep
            write(*,*) wr(i), axr(i), axi(i)
         enddo
         do i = 1, numstep
            write(*,*) ax(i)
         enddo
      else if (dipxyz(m) == 2) then
         write(*,*) 'ay'
         do i = 1, numstep
            write(*,*) wr(i), ayr(i), ayi(i)
         enddo
         do i = 1, numstep
            write(*,*) ay(i)
         enddo
      else if (dipxyz(m) == 3) then
         write(*,*) 'az'
         do i = 1, numstep
            write(*,*) wr(i), azr(i), azi(i)
         enddo
         do i = 1, numstep
            write(*,*) az(i)
         enddo
      endif
   enddo
endif

deallocate(wr)
deallocate(ax,ay,az)
deallocate(intfit)

deallocate(axr,axi,ayr,ayi,azr,azi)

!CLOSE(200)


!t4 = MPI_wtime()
!dt2 = t4 -t3

!write(*,*) myrow, mycol, 'tempo_totale', dt2

END SUBROUTINE
