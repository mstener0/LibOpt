SUBROUTINE Calcalp(context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,&
          &dipxyz,atomtype,numatom,atominv,atomfit,totmaxfit,   &
          &invfit,maxfit,intervals,Ein,Efin,numstep,wi,YN,ANYN,delta,DV)
!numpoint, &
!          &enpoint)


use mpi

implicit none

!input
integer, intent(in) :: context,nprow,npcol,myrow,mycol, total, point
integer, intent(in) :: kdim(point+1), xyz, dipxyz(3)
integer, intent(in) :: atomtype,numatom,atominv(numatom)
integer, intent(in) :: atomfit(numatom),totmaxfit, maxfit(atomtype)
integer, intent(in) :: invfit(atomtype,totmaxfit)
real*8, intent(in) :: intervals(point+1), delta
real*8, intent(inout) :: Ein, Efin, wi!, enpoint(numpoint)
integer, intent(inout) :: numstep
!integer, intent(in) :: numpoint
character*1, intent(in) :: YN, DV, ANYN

integer, parameter :: double = selected_real_kind(14)
real(double), parameter :: eVha = 27.2113961
real(double), parameter :: rotau2cgs = 235.726327d0
integer :: i, j, m, db, ien, k, istat
integer :: iin, ifi

! MPI
integer :: ierr, fa, fh, fs, fd, fmg 
integer :: status(MPI_STATUS_SIZE)
! BLACS/SCALAPACK
integer :: info
external :: blacs_exit, blacs_gridexit, blacs_gridinfo, blacs_get, &
            &blacs_gridinit, blacs_pinfo
external :: descinit, numroc, pdelset, pdlaprnt, indxl2g, pdgemm
external :: pzelset, pzlaprnt, pdgeadd, pzgesv, pzgemv
integer :: numroc, indxl2g

! Matrix
integer :: Mloc, Nloc, NAkloc, ecc, NVloc, Ndiploc, Mdiploc, Ndloc

real(double), allocatable :: Ak(:,:), Dk(:,:), Qr(:,:), Qi(:,:), L(:,:), work(:)
real(double), allocatable :: Mr(:,:), Mi(:,:), S(:,:), Vx(:,:), Vy(:,:), Vz(:,:)
real(double), allocatable :: dxr(:), dxi(:), dyr(:), dyi(:), dzr(:), dzi(:)
real(double), allocatable :: gxr(:), gxi(:), gyr(:), gyi(:), gzr(:), gzi(:)
real(double), allocatable :: nx(:), ny(:), nz(:), wr(:)
real(double), allocatable :: axr(:), axi(:), ayr(:), ayi(:), azr(:), azi(:)
real(double), allocatable :: betxr(:),betxi(:),betyr(:),betyi(:),betzr(:),betzi(:)
real(double), allocatable :: bet2xr(:),bet2xi(:),bet2yr(:),bet2yi(:),bet2zr(:),bet2zi(:)
real(double), allocatable :: bxr(:), bxi(:), byr(:), byi(:), bzr(:), bzi(:)
real(double), allocatable :: vectr(:), vecti(:), dip(:,:), vecAr(:), vecAi(:)
real(double), allocatable :: intfit(:), testr(:), testi(:)
real(double), allocatable :: testran(:), testian(:)
real(double), allocatable :: vectmr(:), vectmi(:)
real(double), allocatable :: qxr(:), qxi(:), qyr(:), qyi(:), qzr(:), qzi(:)

real(double), allocatable :: mag(:,:), Vmx(:,:), Vmy(:,:), Vmz(:,:)


integer, allocatable :: ipvt(:)
integer :: descAk(9), descMat(9), descV(9)
integer :: descVec(9), descdip(9), descvect(9)
integer :: descd(9)
integer :: descsc(9), descscb(9)

real(double) :: aval, val, skr, ski, valr, vali, Ek
real(double) :: tkr, tki, valran, valian
real(double) :: step, nn
real(double) :: scalar(1,1)

integer :: totecc

complex(double) :: sk, tk, w, term1, term2, valc, valc2, tempc
complex(double), allocatable :: Cx(:,:), dx(:), dy(:), dz(:)
complex(double), allocatable :: Cy(:,:), Cz(:,:)
complex(double), allocatable :: ax(:), ay(:), az(:)
complex(double), allocatable :: gx(:), gy(:), gz(:)

complex(double) :: valalp, valbet, alphac, betac
complex(double), allocatable :: intfitcompl(:)
complex(double), allocatable :: dxp(:,:), dyp(:,:), dzp(:,:)
complex(double) :: scalarb(1,3) 

integer(KIND=MPI_OFFSET_KIND) :: offset, rig, col, temp, kiniz

real(double) :: t1, t2, dt, alpha, beta, malpha, t3, t4, dt2

! Initialize a default BLACS context and the processes grid

t3 = MPI_wtime()

alpha = 1.0d0
malpha = -1.0d0
beta = 0.0d0

alphac = CMPLX(alpha,beta)
betac = CMPLX(beta,beta)

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

call descinit(descsc, 1, 1, 1, 1, 0, 0, context, 1, info)
call descinit(descscb, 1, 3, 1, 1, 0, 0, context, 1, info)




if (point < npcol) then
   call descinit(descV, total, point, total/nprow, 1,&
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
dip = 0.0d0

allocate ( vecAr ( Mdiploc ) )
allocate ( vecAi ( Mdiploc ) )
allocate ( testr ( Mdiploc ) )
allocate ( testi ( Mdiploc ) )
allocate ( testran ( Mdiploc ) )
allocate ( testian ( Mdiploc ) )

testr = 0.0d0
testi = 0.0d0
testran = 0.0d0
testian = 0.0d0

if ( totecc < nprow) then
   call descinit(descvect, totecc, 1, 1, 1,&
                & 0, 0, context, max(1,Mdiploc), info)
else
   call descinit(descvect, totecc, 1, totecc/nprow, 1,&
                & 0, 0, context, max(1,Mdiploc), info)
endif

if (3 < npcol) then
   if ( totecc < nprow) then
      call descinit(descdip, totecc, 3, 1, 1,&
                   & 0, 0, context, max(1,Mdiploc), info)
   else
      call descinit(descdip, totecc, 3, totecc/nprow, 1,&
                   & 0, 0, context, max(1,Mdiploc), info)
   endif
else
   if ( totecc < nprow) then
      call descinit(descdip, totecc, 3, 1, 3/npcol,&
                   & 0, 0, context, max(1,Mdiploc), info)
   else
      call descinit(descdip, totecc, 3, totecc/nprow, 3/npcol,&
                   & 0, 0, context, max(1,Mdiploc), info)
   endif
endif


if (DV .eq. 'Y') then
   call MPI_FILE_OPEN(MPI_COMM_WORLD,'dipvel.out',MPI_MODE_RDONLY&
                     &,MPI_INFO_NULL,fd,ierr)

   CALL CalcVxyz(fa,fd,context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,dipxyz,&
             & descV, Mloc, NVloc, Vx, Vy, Vz, Ndiploc, Mdiploc, descdip, dip)

   call MPI_FILE_CLOSE(fd,ierr)
else
   call MPI_FILE_OPEN(MPI_COMM_WORLD,'dipole.out',MPI_MODE_RDONLY&
                     &,MPI_INFO_NULL,fd,ierr)

   CALL CalcVxyz(fa,fd,context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,dipxyz,&
             & descV, Mloc, NVloc, Vx, Vy, Vz, Ndiploc, Mdiploc, descdip, dip)

   call MPI_FILE_CLOSE(fd,ierr)
endif


! magnetic part

allocate ( mag ( Mdiploc, Ndiploc ) )
mag = 0.0d0

allocate( Vmx ( Mloc, NVloc ) )
allocate( Vmy ( Mloc, NVloc ) )
allocate( Vmz ( Mloc, NVloc ) )

Vmx = 0.0d0
Vmy = 0.0d0
Vmz = 0.0d0

call MPI_FILE_OPEN(MPI_COMM_WORLD,'magnetic.out',MPI_MODE_RDONLY&
                  &,MPI_INFO_NULL,fmg,ierr)

CALL CalcVxyz(fa,fmg,context,nprow,npcol,myrow,mycol,total,point,kdim,xyz,dipxyz,&
          & descV, Mloc, NVloc, Vmx, Vmy, Vmz, Ndiploc, Mdiploc, descdip, mag)

call MPI_FILE_CLOSE(fmg,ierr)


Nloc = numroc( total, total/nprow, mycol, 0, npcol )

call descinit(descMat, total, total, total/nprow, total/nprow,&
             & 0, 0, context, max(1,Mloc), info)
call descinit(descVec, total, 1, total/nprow, 1, 0, 0, context&
             &, max(1,Mloc), info)

if (3 < npcol) then
   Ndloc = 1
else
   Ndloc = numroc( 3, 3/npcol, mycol, 0, npcol )
endif

if (3 < npcol) then
   if ( total < nprow) then
      call descinit(descd, total, 3, 1, 1,&
                   & 0, 0, context, max(1,Mloc), info)
   else
      call descinit(descd, total, 3, total/nprow, 1,&
                   & 0, 0, context, max(1,Mloc), info)
   endif
else
   if ( total < nprow) then
      call descinit(descd, total, 3, 1, 3/npcol,&
                   & 0, 0, context, max(1,Mloc), info)
   else
      call descinit(descd, total, 3, total/nprow, 3/npcol,&
                   & 0, 0, context, max(1,Mloc), info)
   endif
endif

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
allocate ( betxr ( numstep ) )
allocate ( betxi ( numstep ) )
allocate ( betyr ( numstep ) )
allocate ( betyi ( numstep ) )
allocate ( betzr ( numstep ) )
allocate ( betzi ( numstep ) )
allocate ( bet2xr ( numstep ) )
allocate ( bet2xi ( numstep ) )
allocate ( bet2yr ( numstep ) )
allocate ( bet2yi ( numstep ) )
allocate ( bet2zr ( numstep ) )
allocate ( bet2zi ( numstep ) )



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

betxr = 0.0d0
betxi = 0.0d0
betyr = 0.0d0
betyi = 0.0d0
betzr = 0.0d0
betzi = 0.0d0
bet2xr = 0.0d0
bet2xi = 0.0d0
bet2yr = 0.0d0
bet2yi = 0.0d0
bet2zr = 0.0d0
bet2zi = 0.0d0

DO ien = 1, numstep
      
   allocate( Qr ( Mloc, Nloc ) )
   allocate( Qi ( Mloc, Nloc ) )
   allocate( Dk ( Mloc, Nloc ) )
   allocate( work( total ) )
   allocate( dxr ( Mloc ) )
   allocate( dxi ( Mloc ) )
   allocate( dyr ( Mloc ) )
   allocate( dyi ( Mloc ) )
   allocate( dzr ( Mloc ) )
   allocate( dzi ( Mloc ) )
   allocate( gxr ( Mloc ) )
   allocate( gxi ( Mloc ) )
   allocate( gyr ( Mloc ) )
   allocate( gyi ( Mloc ) )
   allocate( gzr ( Mloc ) )
   allocate( gzi ( Mloc ) )


   Qr = 0.0d0
   Qi = 0.0d0
   dxr = 0.0d0
   dxi = 0.0d0
   dyr = 0.0d0
   dyi = 0.0d0
   dzr = 0.0d0
   dzi = 0.0d0
   gxr = 0.0d0
   gxi = 0.0d0
   gyr = 0.0d0
   gyi = 0.0d0
   gzr = 0.0d0
   gzi = 0.0d0

   w = CMPLX(wr(ien),wi)

   iin= 1 
   ifi= point

   DO i = 1, point
      if (intervals(i) >= wr(ien)+delta) then
         ifi = i
         exit
      endif
   ENDDO

   kiniz = 0

   DO i = iin, ifi

      ecc = kdim(i+1) - kdim(i)
   
      if (ecc <= 0) CYCLE
      Ek = (intervals(i) + intervals(i+1))/2
      tempc = Ek
      term1 = alphac/(tempc-w)
      term2 = alphac/(tempc+w)

      sk = 2*(term1+term2)
      skr = REAL(sk)
      ski = AIMAG(sk)

      tk = 2*(term1-term2)
      tkr = REAL(tk)
      tki = AIMAG(tk)


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

      CALL pdgemm('N','T',total,total,ecc,alpha,Ak,1,1,descAk,Ak,1,1,descAk,&
                 & beta,Dk,1,1,descMat)

      deallocate ( Ak )

      call pdgeadd('N', total,total,skr,Dk,1,1,descMat,alpha,Qr,1,1,descMat)
      call pdgeadd('N', total,total,ski,Dk,1,1,descMat,alpha,Qi,1,1,descMat)

      call pdgeadd('N', total,1,skr,Vx,1,i,descV,alpha,dxr,1,1,descVec)
      call pdgeadd('N', total,1,ski,Vx,1,i,descV,alpha,dxi,1,1,descVec)
      call pdgeadd('N', total,1,skr,Vy,1,i,descV,alpha,dyr,1,1,descVec)
      call pdgeadd('N', total,1,ski,Vy,1,i,descV,alpha,dyi,1,1,descVec)
      call pdgeadd('N', total,1,skr,Vz,1,i,descV,alpha,dzr,1,1,descVec)
      call pdgeadd('N', total,1,ski,Vz,1,i,descV,alpha,dzi,1,1,descVec)

      call pdgeadd('N', total,1,-tkr,Vmx,1,i,descV,alpha,gxi,1,1,descVec)
      call pdgeadd('N', total,1,tki,Vmx,1,i,descV,alpha,gxr,1,1,descVec)
      call pdgeadd('N', total,1,-tkr,Vmy,1,i,descV,alpha,gyi,1,1,descVec)
      call pdgeadd('N', total,1,tki,Vmy,1,i,descV,alpha,gyr,1,1,descVec)
      call pdgeadd('N', total,1,-tkr,Vmz,1,i,descV,alpha,gzi,1,1,descVec)
      call pdgeadd('N', total,1,tki,Vmz,1,i,descV,alpha,gzr,1,1,descVec)

   enddo

   allocate ( dx ( Mloc ) )
   allocate ( dy ( Mloc ) )
   allocate ( dz ( Mloc ) )

   dx = CMPLX (dxr, dxi )
   dy = CMPLX (dyr, dyi )
   dz = CMPLX (dzr, dzi )

   allocate ( gx ( Mloc ) )
   allocate ( gy ( Mloc ) )
   allocate ( gz ( Mloc ) )

   gx = CMPLX (gxr, gxi )
   gy = CMPLX (gyr, gyi )
   gz = CMPLX (gzr, gzi )

! MODIFICHE
! -------------------------------------------------------------------------


   allocate (dxp ( Mloc, Ndloc) )
   allocate (dyp ( Mloc, Ndloc) )
   allocate (dzp ( Mloc, Ndloc) )

   dxp = CMPLX(0.0d0,0.0d0)
   dyp = CMPLX(0.0d0,0.0d0)
   dzp = CMPLX(0.0d0,0.0d0)

   call pzgeadd('N', total,1,alphac,dx,1,1,descVec,betac,dxp,1,1,descd)
   call pzgeadd('N', total,1,alphac,dy,1,1,descVec,betac,dyp,1,1,descd)
   call pzgeadd('N', total,1,alphac,dz,1,1,descVec,betac,dzp,1,1,descd)

   call pzgeadd('N', total,1,alphac,gx,1,1,descVec,betac,dxp,1,3,descd)
   call pzgeadd('N', total,1,alphac,gy,1,1,descVec,betac,dyp,1,3,descd)
   call pzgeadd('N', total,1,alphac,gz,1,1,descVec,betac,dzp,1,3,descd)

   if (YN == 'Y') then

      call pzgeadd('N', total,1,alpha,intfitcompl,1,1,descVec,beta,dxp,1,2,descd)
      call pzgeadd('N', total,1,alpha,intfitcompl,1,1,descVec,beta,dyp,1,2,descd)
      call pzgeadd('N', total,1,alpha,intfitcompl,1,1,descVec,beta,dzp,1,2,descd)

   endif

   deallocate(dxr,dxi,dyr,dyi,dzr,dzi)
   deallocate(gxr,gxi,gyr,gyi,gzr,gzi)
   deallocate(Dk)

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

         rig = INDXL2G(descMat(5)+1, descMat(5), myrow, 0, nprow)
         offset = (temp + rig -1) * db
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

   allocate( Mr ( Mloc, Nloc ) )
   allocate( Mi ( Mloc, Nloc ) )

   CALL pdgemm('N','N',total,total,total,alpha,Qr,1,1,descMat,L,1,1,descMat,&
              & beta,Mr,1,1,descMat)

   CALL pdgemm('N','N',total,total,total,alpha,Qi,1,1,descMat,L,1,1,descMat,&
              & beta,Mi,1,1,descMat)


   deallocate(Qr, Qi)

   allocate( S ( Mloc, Nloc ) )

   call MPI_FILE_OPEN(MPI_COMM_WORLD,'matrixS.out',MPI_MODE_RDONLY&
                     &,MPI_INFO_NULL,fs,ierr)

   DO i = 1, Nloc
      col = INDXL2G(i, total/nprow, mycol, 0, npcol)
      temp = (col - 1) * total
      IF (descMat(9)>descMat(5)) then
         rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
         offset = (temp + rig -1) * db
         CALL MPI_FILE_READ_AT(fs, offset, S(1,i), descMat(5),&
                              & MPI_DOUBLE_PRECISION,status, ierr)
         rig = INDXL2G(descMat(5)+1, descMat(5), myrow, 0, nprow)
         offset = (temp + rig -1) * db

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

   allocate( Cx ( Mloc, Nloc ) )
   allocate( Cy ( Mloc, Nloc ) )
   allocate( Cz ( Mloc, Nloc ) )

   do j = 1, Nloc
      do m = 1, Mloc
         valr = S(m,j)   +Mr(m,j)       
         vali = 0.0d0    +Mi(m,j)      
         Cx(m,j) = CMPLX ( valr , vali )
      enddo
   enddo
  
   Cy = Cx
   Cz = Cy

   deallocate(Mr,Mi,S)

   t1 = MPI_wtime()

   do m = 1, xyz

testran = 0.0d0
testian = 0.0d0



      if (dipxyz(m) == 1) then


         allocate (bxr ( Mloc ) )
         allocate (bxi ( Mloc ) )

         allocate (qxr ( Mloc ) )
         allocate (qxi ( Mloc ) )


         call pzgesv(total,3,Cx,1,1,descMat,ipvt,dxp,1,1,descd,info)

         if (YN == 'N') then

!            call pzgesv(total,3,Cx,1,1,descMat,ipvt,dxp,1,1,descd,info)

            call pzgeadd('N',total,1,alphac,dxp,1,1,descd,betac,dx,1,1,descVec)

            bxr = REAL(dx)
            bxi = AIMAG(dx)

!-----------------------------------------------------------------------------------
!DICROISMO II

            call pzgeadd('N',total,1,alphac,dxp,1,3,descd,betac,gx,1,1,descVec)

            qxr = REAL(gx)
            qxi = AIMAG(gx)

!----------------------------------------------------------------------------------

         else

!            call pzgesv(total,3,Cx,1,1,descMat,ipvt,dxp,1,1,descd,info)

            valalp = betac

            valbet = betac

            call pzgemv('T',total,3,alphac,dxp,1,1,descd,intfitcompl,1,1,descVec,&
                    & 1,betac,scalarb,1,1,descscb,1)

            call pzelget('A', 'D', valalp, scalarb, 1, 1, descscb)
            call pzelget('A', 'D', valbet, scalarb, 1, 2, descscb)


            valc = valalp/valbet

            call pzgeadd('N',total,1,-valc,dxp,1,2,descd,alpha,dxp,1,1,descd)

            call pzgeadd('N',total,1,alphac,dxp,1,1,descd,betac,dx,1,1,descVec)

            bxr = REAL(dx)
            bxi = AIMAG(dx)
            
!----------------------------------------------------------------------------------
!DICROISMO II

            valalp = betac

            call pzelget('A', 'D', valalp, scalarb, 1, 3, descscb)

            valc = valalp/valbet

            call pzgeadd('N',total,1,-valc,dxp,1,2,descd,alpha,dxp,1,3,descd)
            call pzgeadd('N', total,1,alphac,dxp,1,3,descd,betac,gx,1,1,descVec)

            qxr = REAL(gx)
            qxi = AIMAG(gx)

!----------------------------------------------------------------------------------

         endif

         if (ANYN == 'Y') then
            OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
                 & POSITION='APPEND',IOSTAT=istat)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'point', ien, 'energy', wr(ien)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'bx'
            do k = 1, total
               valran = beta
               valian = beta
               call pdelget('A', 'D', valran, bxr, k, 1, descvect)
               call pdelget('A', 'D', valian, bxi, k, 1, descvect)
               if (myrow == 0 .and. mycol ==0) then
                  write(200,*) valran, valian
               endif
            enddo
            CLOSE(200)
         endif

!         do i = 1, numpoint
!            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
!            if (ANYN == 'Y') then
!               do j = 0, nprow-1
!                  if (myrow == j .and. mycol ==0) then
!                     OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!                     & POSITION='APPEND',IOSTAT=istat)
!                     if (j == 0) write(200,*) 'point', ien, 'energy', wr(ien) !enpoint(i)
!                     if (j == 0) write(200,*) 'bx'
!                     do k = 1, Mloc
!                        write(200,*) bxr(k), bxi(k)
!                     enddo
!                     CLOSE(200)
!                  endif
!                  call MPI_Barrier( mpi_comm_world, ierr)
!               enddo
!            endif
!            endif
!         enddo

         allocate (vectr ( Mloc ) )
         allocate (vecti ( Mloc ) )

         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bxr,1,1,descVec,&
                    & 1,beta,vectr,1,1,descVec,1)
         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bxi,1,1,descVec,&
                    & 1,beta,vecti,1,1,descVec,1)

!-------------------------------------------------------------------------------------
! DICROISMO II

         allocate (vectmr ( Mloc ) )
         allocate (vectmi ( Mloc ) )

         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,qxr,1,1,descVec,&
                    & 1,beta,vectmr,1,1,descVec,1)
         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,qxi,1,1,descVec,&
                    & 1,beta,vectmi,1,1,descVec,1)

!-------------------------------------------------------------------------------------

         kiniz = 0

         DO i = iin, ifi

            ecc = kdim(i+1) - kdim(i)
            if (ecc <= 0) CYCLE

            Ek = (intervals(i) + intervals(i+1))/2
            tempc = Ek
            term1 = alphac/(tempc-w)
            term2 = alphac/(tempc+w)
            sk = 2*(term1+term2)

            skr = REAL(sk)
            ski = AIMAG(sk)

            tk = 2*(term1-term2)
            tkr = REAL(tk)
            tki = AIMAG(tk)

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
            kiniz = kiniz + ecc*total


            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectr,1,1,descVec,&
                       & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vecti,1,1,descVec,&
                       & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

!            deallocate ( Ak )

            call pdgeadd('N',ecc,1,alpha,dip,kdim(i)+1,1,descdip,malpha,vecAr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,skr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-skr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            valr = beta  

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,1,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)
     
            vali = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,1,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            axr(ien) = axr(ien) + valr      
            axi(ien) = axi(ien) + vali 

            call pdgeadd('N',ecc,1,alpha,testr,kdim(i)+1,1,descvect,beta,testran,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,alpha,testi,kdim(i)+1,1,descvect,beta,testian,kdim(i)+1,1,descvect)


!----------------------------------------------------------------------------------------------------------
! DICROISMO

            call pdgeadd('N',ecc,1,tkr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,tki,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-tkr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,tki,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            vali = beta

            CALL pdgemv('T',ecc,1,-alpha,mag,kdim(i)+1,1,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,mag,kdim(i)+1,1,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            betxr(ien) = betxr(ien) + valr
            betxi(ien) = betxi(ien) + vali

!-------------------------------------------------------------------------------------------------------
!DICROISMO II

            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectmr,1,1,descVec,&
                       & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectmi,1,1,descVec,&
                       & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

            deallocate ( Ak )            

            call pdgeadd('N',ecc,1,-skr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-skr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-ski,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            call pdgeadd('N',ecc,1,tki,mag,kdim(i)+1,1,descdip,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-tkr,mag,kdim(i)+1,1,descdip,alpha,testi,kdim(i)+1,1,descvect)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,1,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            vali = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,1,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            bet2xr(ien) = bet2xr(ien) + valr
            bet2xi(ien) = bet2xi(ien) + vali

!-------------------------------------------------------------------------------------------------------

         enddo

         if (ANYN == 'Y') then
            OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
                 & POSITION='APPEND',IOSTAT=istat)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'Rirx'
            do k = 1, totecc
               valran = beta
               valian = beta
               call pdelget('A', 'D', valran, testran, k, 1, descvect)
               call pdelget('A', 'D', valian, testian, k, 1, descvect)
               if (myrow == 0 .and. mycol ==0) then
                  write(200,*) valran, valian
               endif                  
            enddo
            CLOSE(200)
         endif




!         do i = 1, numpoint
!            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
!            if (ANYN == 'Y') then
!               do j = 0, nprow-1 
!                  if (myrow == j .and. mycol ==0) then
!                     OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!                          & POSITION='APPEND',IOSTAT=istat)
!                     if (j == 0) write(200,*) 'Rirx'
!                     do k = 1, Mdiploc
!                        write(200,*) testran(k), testian(k)
!                     enddo
!                     CLOSE(200)
!                  endif
!                  call MPI_Barrier( mpi_comm_world, ierr)
!               enddo
!            endif
!            endif
!         enddo
 
         deallocate (vectr,vecti)
         deallocate (vectmr,vectmi)

         valr = beta

         CALL pdgemv('T',total,1,alpha,bxr,1,1,descVec,nx,1,1,descVec,&
                    & 1,beta,scalar,1,1,descsc,1)

         call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

         vali = beta

         CALL pdgemv('T',total,1,alpha,bxi,1,1,descVec,nx,1,1,descVec,&
                    & 1,beta,scalar,1,1,descsc,1)

         call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

         deallocate(bxr,bxi,qxr,qxi)

         ax(ien) = ax(ien) + CMPLX(valr, vali)

      else if (dipxyz(m) == 2 ) then
 
         allocate (byr ( Mloc ) )
         allocate (byi ( Mloc ) )

         allocate (qyr ( Mloc ) )
         allocate (qyi ( Mloc ) )


         call pzgesv(total,3,Cy,1,1,descMat,ipvt,dyp,1,1,descd,info)

         if (YN == 'N' ) then

!            call pzgesv(total,1,Cy,1,1,descMat,ipvt,dyp,1,1,descd,info)

            call pzgeadd('N',total,1,alphac,dyp,1,1,descd,betac,dy,1,1,descVec)

            byr = REAL(dy)
            byi = AIMAG(dy)

!-----------------------------------------------------------------------------------
!DICROISMO II

            call pzgeadd('N',total,1,alphac,dyp,1,3,descd,betac,gy,1,1,descVec)

            qyr = REAL(gy)
            qyi = AIMAG(gy)

!----------------------------------------------------------------------------------

         else

!            call pzgesv(total,2,Cy,1,1,descMat,ipvt,dyp,1,1,descd,info)

            valalp = betac

            valbet = betac
 
            call pzgemv('T',total,3,alphac,dyp,1,1,descd,intfitcompl,1,1,descVec,&
                    & 1,betac,scalarb,1,1,descscb,1)

            call pzelget('A', 'D', valalp, scalarb, 1, 1, descscb)
            call pzelget('A', 'D', valbet, scalarb, 1, 2, descscb)

            valc = valalp/valbet

            call pzgeadd('N',total,1,-valc,dyp,1,2,descd,alpha,dyp,1,1,descd)

            call pzgeadd('N',total,1,alphac,dyp,1,1,descd,betac,dy,1,1,descVec)

            byr = REAL(dy)
            byi = AIMAG(dy)

!----------------------------------------------------------------------------------
!DICROISMO II

            call pzelget('A', 'D', valalp, scalarb, 1, 3, descscb)

            valc = valalp/valbet

            call pzgeadd('N',total,1,-valc,dyp,1,2,descd,alpha,dyp,1,3,descd)
            call pzgeadd('N', total,1,alphac,dyp,1,3,descd,betac,gy,1,1,descVec)

            qyr = REAL(gy)
            qyi = AIMAG(gy)

!----------------------------------------------------------------------------------

         endif

         if (ANYN == 'Y') then
            OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
                 & POSITION='APPEND',IOSTAT=istat)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'point', ien, 'energy', wr(ien)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'by'
            do k = 1, total
               valran = beta
               valian = beta
               call pdelget('A', 'D', valran, byr, k, 1, descvect)
               call pdelget('A', 'D', valian, byi, k, 1, descvect)
               if (myrow == 0 .and. mycol ==0) then
                  write(200,*) valran, valian
               endif
            enddo
            CLOSE(200)
         endif

!         do i = 1, numpoint
!            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
!            if (ANYN == 'Y') then
!               do j = 0, nprow-1
!                  if (myrow == j .and. mycol ==0) then
!                     OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!                     & POSITION='APPEND',IOSTAT=istat)
!                     if (j == 0) write(200,*) 'point', ien, 'energy', wr(ien) !enpoint(i)
!                     if (j == 0) write(200,*) 'by'
!                     do k = 1, Mloc
!                        write(200,*) byr(k), byi(k)
!                     enddo
!                     CLOSE(200)
!                  endif
!                  call MPI_Barrier( mpi_comm_world, ierr)
!               enddo
!            endif
!            endif
!         enddo


         allocate (vectr ( Mloc ) )
         allocate (vecti ( Mloc ) )

         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,byr,1,1,descVec,&
                    & 1,beta,vectr,1,1,descVec,1)
         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,byi,1,1,descVec,&
                    & 1,beta,vecti,1,1,descVec,1)

!-------------------------------------------------------------------------------------
! DICROISMO II

         allocate (vectmr ( Mloc ) )
         allocate (vectmi ( Mloc ) )

         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,qyr,1,1,descVec,&
                    & 1,beta,vectmr,1,1,descVec,1)
         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,qyi,1,1,descVec,&
                    & 1,beta,vectmi,1,1,descVec,1)

!-------------------------------------------------------------------------------------

         kiniz = 0

         DO i = iin, ifi

            ecc = kdim(i+1) - kdim(i)
            if (ecc <= 0) CYCLE

            Ek = (intervals(i) + intervals(i+1))/2
            tempc = Ek
            term1 = alphac/(tempc-w)
            term2 = alphac/(tempc+w)
            sk = 2*(term1+term2)


            skr = REAL(sk)
            ski = AIMAG(sk)

            tk = 2*(term1-term2)
            tkr = REAL(tk)
            tki = AIMAG(tk)

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
            kiniz = kiniz + ecc*total


            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectr,1,1,descVec,&
                       & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vecti,1,1,descVec,&
                       & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

!            deallocate ( Ak )

            call pdgeadd('N',ecc,1,alpha,dip,kdim(i)+1,2,descdip,malpha,vecAr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,skr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-skr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            valr = beta   

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,2,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            vali = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,2,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            ayr(ien) = ayr(ien) + valr      
            ayi(ien) = ayi(ien) + vali   

            call pdgeadd('N',ecc,1,alpha,testr,kdim(i)+1,1,descvect,beta,testran,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,alpha,testi,kdim(i)+1,1,descvect,beta,testian,kdim(i)+1,1,descvect)

!----------------------------------------------------------------------------------------------------------
! DICROISMO

            call pdgeadd('N',ecc,1,tkr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,tki,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-tkr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,tki,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            vali = beta

            CALL pdgemv('T',ecc,1,-alpha,mag,kdim(i)+1,2,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,mag,kdim(i)+1,2,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            betyr(ien) = betyr(ien) + valr
            betyi(ien) = betyi(ien) + vali

!-------------------------------------------------------------------------------------------------------
!DICROISMO II

            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectmr,1,1,descVec,&
                       & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectmi,1,1,descVec,&
                       & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

            deallocate ( Ak )

            call pdgeadd('N',ecc,1,-skr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-skr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-ski,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            call pdgeadd('N',ecc,1,tki,mag,kdim(i)+1,2,descdip,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-tkr,mag,kdim(i)+1,2,descdip,alpha,testi,kdim(i)+1,1,descvect)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,2,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            vali = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,2,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            bet2yr(ien) = bet2yr(ien) + valr
            bet2yi(ien) = bet2yi(ien) + vali

!-------------------------------------------------------------------------------------------------------

         enddo

         if (ANYN == 'Y') then
            OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
                 & POSITION='APPEND',IOSTAT=istat)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'Riry'
            do k = 1, totecc
               valran = beta
               valian = beta
               call pdelget('A', 'D', valran, testran, k, 1, descvect)
               call pdelget('A', 'D', valian, testian, k, 1, descvect)
               if (myrow == 0 .and. mycol ==0) then
                  write(200,*) valran, valian
               endif   
            enddo
            CLOSE(200)
         endif



!         do i = 1, numpoint
!            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
!            if (ANYN == 'Y') then
!               do j = 0, nprow-1
!                  if (myrow == j .and. mycol ==0) then
!                     OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!                          & POSITION='APPEND',IOSTAT=istat)
!                     if (j == 0) write(200,*) 'Riry'
!                     do k = 1, Mdiploc
!                        write(200,*) testran(k), testian(k)
!                     enddo
!                     CLOSE(200)
!                  endif
!                  call MPI_Barrier( mpi_comm_world, ierr)
!               enddo
!            endif
!            endif
!         enddo


         deallocate (vectr,vecti)

         deallocate (vectmr,vectmi)

         valr = beta

         CALL pdgemv('T',total,1,alpha,byr,1,1,descVec,ny,1,1,descVec,&
                    & 1,beta,scalar,1,1,descsc,1)

         call pdelget('A', 'D', valr, scalar, 1, 1, descsc)


         vali = beta

         CALL pdgemv('T',total,1,alpha,byi,1,1,descVec,ny,1,1,descVec,&
                    & 1,beta,scalar,1,1,descsc,1)

         call pdelget('A', 'D', vali, scalar, 1, 1, descsc)


         deallocate(byr,byi,qyr,qyi)

         ay(ien) = ay(ien) + CMPLX(valr, vali)

      else if (dipxyz(m) == 3 ) then

         allocate (bzr ( Mloc ) )
         allocate (bzi ( Mloc ) )

         allocate (qzr ( Mloc ) )
         allocate (qzi ( Mloc ) )

         call pzgesv(total,3,Cz,1,1,descMat,ipvt,dzp,1,1,descd,info)

         if (YN == 'N') then

!            call pzgesv(total,1,Cz,1,1,descMat,ipvt,dzp,1,1,descd,info)

            call pzgeadd('N',total,1,alphac,dzp,1,1,descd,betac,dz,1,1,descVec)

            bzr = REAL(dz)
            bzi = AIMAG(dz)

!-----------------------------------------------------------------------------------
!DICROISMO II

            call pzgeadd('N',total,1,alphac,dzp,1,3,descd,betac,gz,1,1,descVec)

            qzr = REAL(gz)
            qzi = AIMAG(gz)

!----------------------------------------------------------------------------------

         else

!            call pzgesv(total,2,Cz,1,1,descMat,ipvt,dzp,1,1,descd,info)

            valalp = betac
            valbet = betac

            call pzgemv('T',total,3,alphac,dzp,1,1,descd,intfitcompl,1,1,descVec,&
                    & 1,betac,scalarb,1,1,descscb,1)

            call pzelget('A', 'D', valalp, scalarb, 1, 1, descscb)
            call pzelget('A', 'D', valbet, scalarb, 1, 2, descscb)

            valc = valalp/valbet

            call pzgeadd('N',total,1,-valc,dzp,1,2,descd,alpha,dzp,1,1,descd)
            call pzgeadd('N', total,1,alphac,dzp,1,1,descd,betac,dz,1,1,descVec)

            bzr = REAL(dz)
            bzi = AIMAG(dz)


!----------------------------------------------------------------------------------
!DICROISMO II

            call pzelget('A', 'D', valalp, scalarb, 1, 3, descscb)
 
            valc = valalp/valbet

            call pzgeadd('N',total,1,-valc,dzp,1,2,descd,alpha,dzp,1,3,descd)
            call pzgeadd('N', total,1,alphac,dzp,1,3,descd,betac,gz,1,1,descVec)

            qzr = REAL(gz)
            qzi = AIMAG(gz)

!----------------------------------------------------------------------------------

         endif

         if (ANYN == 'Y') then
            OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
                 & POSITION='APPEND',IOSTAT=istat)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'point', ien, 'energy', wr(ien)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'bz'
            do k = 1, total
               valran = beta
               valian = beta
               call pdelget('A', 'D', valran, bzr, k, 1, descvect)
               call pdelget('A', 'D', valian, bzi, k, 1, descvect)
               if (myrow == 0 .and. mycol ==0) then
                  write(200,*) valran, valian
               endif
            enddo
            CLOSE(200)
         endif


!         do i = 1, numpoint
!            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
!            if (ANYN == 'Y') then
!               do j = 0, nprow-1
!                  if (myrow == j .and. mycol ==0) then
!                     OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!                          & POSITION='APPEND',IOSTAT=istat)
!                     if (j == 0) write(200,*) 'point', ien, 'energy', wr(ien)!enpoint(i)
!                     if (j == 0) write(200,*) 'bz'
!                     do k = 1, Mloc
!                        write(200,*) bzr(k), bzi(k)
!                     enddo
!                     CLOSE(200)
!                  endif
!                  call MPI_Barrier( mpi_comm_world, ierr)
!               enddo
!            endif
!            endif
!         enddo

         allocate (vectr ( Mloc ) )
         allocate (vecti ( Mloc ) )         

         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bzr,1,1,descVec,&
                    & 1,beta,vectr,1,1,descVec,1)
         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,bzi,1,1,descVec,&
                    & 1,beta,vecti,1,1,descVec,1) 

!-------------------------------------------------------------------------------------
! DICROISMO II

         allocate (vectmr ( Mloc ) )
         allocate (vectmi ( Mloc ) )

         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,qzr,1,1,descVec,&
                    & 1,beta,vectmr,1,1,descVec,1)
         CALL pdgemv('N',total,total,alpha,L,1,1,descMat,qzi,1,1,descVec,&
                    & 1,beta,vectmi,1,1,descVec,1)

!-------------------------------------------------------------------------------------

         kiniz = 0

         DO i = iin, ifi
   
            ecc = kdim(i+1) - kdim(i)   
            if (ecc <= 0) CYCLE

            Ek = (intervals(i) + intervals(i+1))/2
            tempc = Ek
            term1 = alphac/(tempc-w)
            term2 = alphac/(tempc+w)
            sk = 2*(term1+term2)

            skr = REAL(sk)
            ski = AIMAG(sk)

            tk = 2*(term1-term2)
            tkr = REAL(tk)
            tki = AIMAG(tk)

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
                  CALL MPI_FILE_READ_AT(fa,offset,Ak(1+descAk(5),j),&
                                       &descAk(9)-descAk(5), MPI_DOUBLE_PRECISION,&
                                       &status, ierr)
               ELSE
                  rig = INDXL2G(1, total/nprow, myrow, 0, nprow)
                  offset = (temp+rig-1)*db
                  CALL MPI_FILE_READ_AT(fa,offset,Ak(1,j),descAk(5),&
                                       & MPI_DOUBLE_PRECISION,status,ierr)
               ENDIF
            ENDDO
            kiniz = kiniz + ecc*total
    
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectr,1,1,descVec,&
                       & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vecti,1,1,descVec,&
                       & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

!            deallocate ( Ak )

            call pdgeadd('N',ecc,1,alpha,dip,kdim(i)+1,3,descdip,malpha,vecAr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,skr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-skr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,3,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            vali = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,3,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            azr(ien) = azr(ien) + valr      
            azi(ien) = azi(ien) + vali       

            call pdgeadd('N',ecc,1,alpha,testr,kdim(i)+1,1,descvect,beta,testran,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,alpha,testi,kdim(i)+1,1,descvect,beta,testian,kdim(i)+1,1,descvect)

!----------------------------------------------------------------------------------------------------------
! DICROISMO I

            call pdgeadd('N',ecc,1,tkr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,tki,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-tkr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,tki,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            vali = beta


            CALL pdgemv('T',ecc,1,-alpha,mag,kdim(i)+1,3,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,mag,kdim(i)+1,3,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            betzr(ien) = betzr(ien) + valr
            betzi(ien) = betzi(ien) + vali

!-------------------------------------------------------------------------------------------------------
!DICROISMO II

            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectmr,1,1,descVec,&
                       & 1,beta,vecAr,kdim(i)+1,1,descvect,1)
            CALL pdgemv('T',total,ecc,alpha,Ak,1,1,descAk,vectmi,1,1,descVec,&
                       & 1,beta,vecAi,kdim(i)+1,1,descvect,1)

            deallocate ( Ak )

            call pdgeadd('N',ecc,1,-skr,vecAr,kdim(i)+1,1,descvect,beta,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,ski,vecAi,kdim(i)+1,1,descvect,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-skr,vecAi,kdim(i)+1,1,descvect,beta,testi,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-ski,vecAr,kdim(i)+1,1,descvect,alpha,testi,kdim(i)+1,1,descvect)

            call pdgeadd('N',ecc,1,tki,mag,kdim(i)+1,3,descdip,alpha,testr,kdim(i)+1,1,descvect)
            call pdgeadd('N',ecc,1,-tkr,mag,kdim(i)+1,3,descdip,alpha,testi,kdim(i)+1,1,descvect)

            valr = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,3,descdip,testr,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', valr, scalar, 1, 1, descsc)

            vali = beta

            CALL pdgemv('T',ecc,1,alpha,dip,kdim(i)+1,3,descdip,testi,kdim(i)+1,1,descvect,&
                       & 1,beta,scalar,1,1,descsc,1)

            call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

            bet2zr(ien) = bet2zr(ien) + valr
            bet2zi(ien) = bet2zi(ien) + vali

!-------------------------------------------------------------------------------------------------------

         enddo

         if (ANYN == 'Y') then
            OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
                 & POSITION='APPEND',IOSTAT=istat)
            if (myrow == 0 .and. mycol ==0) write(200,*) 'Rirz'
            do k = 1, totecc
               valran = beta
               valian = beta
               call pdelget('A', 'D', valran, testran, k, 1, descvect)
               call pdelget('A', 'D', valian, testian, k, 1, descvect)
               if (myrow == 0 .and. mycol ==0) then
                  write(200,*) valran, valian
               endif
            enddo
            CLOSE(200)
         endif



!         do i = 1, numpoint
!            if ( enpoint(i) >= wr(ien) .and. enpoint(i) < wr(ien+1)) then
!            if (ANYN == 'Y') then
!               do j = 0, nprow-1
!                  if (myrow == j .and. mycol ==0) then
!                     OPEN (UNIT=200, FILE= 'analysis', STATUS='OLD',ACTION='WRITE', &
!                          & POSITION='APPEND',IOSTAT=istat)
!                     if (j == 0) write(200,*) 'Rirz'
!                     do k = 1, Mdiploc
!                        write(200,*) testran(k), testian(k)
!                     enddo
!                     CLOSE(200)
!                  endif
!                  call MPI_Barrier( mpi_comm_world, ierr)
!               enddo
!            endif
!            endif
!         enddo

         deallocate (vectr,vecti)

         deallocate (vectmr,vectmi)    

         valr = beta

         CALL pdgemv('T',total,1,alpha,bzr,1,1,descVec,nz,1,1,descVec,&
                    & 1,beta,scalar,1,1,descsc,1)

         call pdelget('A', 'D', valr, scalar, 1, 1, descsc)


         vali = beta

         CALL pdgemv('T',total,1,alpha,bzi,1,1,descVec,nz,1,1,descVec,&
                    & 1,beta,scalar,1,1,descsc,1)

         call pdelget('A', 'D', vali, scalar, 1, 1, descsc)

         deallocate(bzr,bzi,qzr,qzi)

         az(ien) = az(ien) + CMPLX(valr, vali)

      endif
   enddo

   t2 = MPI_wtime()
   dt = t2 -t1

   deallocate(L)
   deallocate(dx,dy,dz)
   deallocate(gx,gy,gz)
   deallocate(dxp,dyp,dzp)

   deallocate( Cx, Cy, Cz )
   deallocate( work )

enddo

deallocate(mag,dip,vecAr,vecAi,testr,testi)
deallocate(testran,testian)

deallocate(nx, ny, nz)
if (YN=='Y') then
   deallocate(intfitcompl)
endif

deallocate(ipvt)

deallocate(Vx,Vy,Vz)
deallocate(Vmx,Vmy,Vmz)

!chiudere ciclo energie

call MPI_FILE_CLOSE(fa,ierr)

if (myrow==0 .and. mycol==0) then
   do m = 1, xyz
      if (dipxyz(m) == 1) then
         write(*,*) 'ax'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, axr(i), axi(i)*wr(i)*(2./3.)*wi
         enddo
         do i = 1, numstep
            write(*,*) REAL(ax(i)), AIMAG(ax(i)), AIMAG(ax(i))*wr(i)*wi
         enddo
      else if (dipxyz(m) == 2) then
         write(*,*) 'ay'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, ayr(i), ayi(i)*wr(i)*(2./3.)*wi
         enddo
         do i = 1, numstep
            write(*,*) REAL(ay(i)), AIMAG(ay(i)), AIMAG(ay(i))*wr(i)*wi
         enddo
      else if (dipxyz(m) == 3) then
         write(*,*) 'az'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, azr(i), azi(i)*wr(i)*(2./3.)*wi
         enddo
         do i = 1, numstep
            write(*,*) REAL(az(i)), AIMAG(az(i)), AIMAG(az(i))*wr(i)*wi
         enddo
      endif
   enddo
endif

!---------------------------------------------------------------------------------------
!DICROISMO


if (myrow==0 .and. mycol==0) then
   do m = 1, xyz
      if (dipxyz(m) == 1) then
         write(*,*) 'betx'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, betxr(i)*wi*rotau2cgs
         enddo
      else if (dipxyz(m) == 2) then
         write(*,*) 'bety'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, betyr(i)*wi*rotau2cgs
         enddo
      else if (dipxyz(m) == 3) then
         write(*,*) 'betz'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, betzr(i)*wi*rotau2cgs
         enddo
      endif
   enddo
endif

!---------------------------------------------------------------------------------------
!DICROISMO II

if (myrow==0 .and. mycol==0) then
   do m = 1, xyz
      if (dipxyz(m) == 1) then
         write(*,*) 'bet2x'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, bet2xr(i)*wi*rotau2cgs
         enddo
      else if (dipxyz(m) == 2) then
         write(*,*) 'bet2y'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, bet2yr(i)*wi*rotau2cgs
         enddo
      else if (dipxyz(m) == 3) then
         write(*,*) 'bet2z'
         do i = 1, numstep
            write(*,3000) wr(i), wr(i)*eVha, bet2zr(i)*wi*rotau2cgs
         enddo
      endif
   enddo
endif

!---------------------------------------------------------------------------------------

deallocate(wr)
deallocate(ax,ay,az)
deallocate(intfit)

deallocate(axr,axi,ayr,ayi,azr,azi)
deallocate(betxr,betxi,betyr,betyi,betzr,betzi)
deallocate(bet2xr,bet2xi,bet2yr,bet2yi,bet2zr,bet2zi)

3000 format (1x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8)

END SUBROUTINE
