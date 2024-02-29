Subroutine CalcAlp()

   use KF
   use PoltddftModule
   use ppCommonsModule   
   use Vartypes
   use ExcitSOVars
   use ExcitOpenMod
   use SymmetryInfo
   use DimensionsModule
   use DistributedMatrixModule
   use ScaMatrixArrayModule
   use SharedArraysModule
   use adf_blacs
   use ADFFilesModule
   use PairsModule
   use MasterModule

   implicit none 

   integer(KINT) :: ien, ifi, i, ki, kd, kiEnd, m, k, n
   integer(KINT) :: iu61, iu62, iu63, iu21
   integer(KINT) :: context
   integer(KINT) :: block
   integer(KINT) :: info

   real(KREAL)   :: skr, ski, tkr, tki

   complex(KREAL) :: sk, tk, w
   complex(KREAL) :: tempc, term1, term2
   complex(KREAL) :: vala, valb, valc, vald

   real(KREAL), pointer :: sharedAk(:,:)

   real(KREAL), allocatable    :: intfit(:)
   complex(KREAL), allocatable :: intfitc(:)
   real(KREAL), allocatable    :: vecAr(:), vecAi(:), testr(:), testi(:)
   type(DistributedMatrixType) :: Qr, Qi, Dk
   type(DistributedMatrixType) :: S, L, Mr, Mi
   type(DistributedMatrixType) :: C
   complex(KREAL), allocatable :: dc(:,:), dp(:,:), dm(:,:) !analysis 2
   real(KREAL), allocatable    :: vectr(:), vecti(:)
   real(KREAL), allocatable    :: br(:), bi(:), gr(:), gi(:)

   character(LCHARS) :: idxlabel, kstr

   complex(KREAL), parameter :: one_c = (1.0_KREAL, 0.0_KREAL)
   complex(KREAL), parameter :: zero_c  = (0.0_KREAL, 0.0_KREAL)

   logical, external :: master

   call timers ('Alp')
   call timers ('pp')

   call GetBLACSContextAndBlocksize(gDims%nsfos, context, block)

   allocate(intfit(gDims%nsfos))
   intfit = 0.0_KREAL

   call kfopfl(iu62, 'TAPE11')
   call kfread(iu62, 'intfit%intfit', intfit)
   call kfclfl(iu62)

   call kfopfl (iu63, 'TAPE63')
   call kfopsc (iu63, 'MatrixAk')
   call kfopfl (iu61, 'TAPE61')
   call kfopsc (iu61, 'Matrix')

   if (master()) then
      call CreateKFPolTDDFT
      call kfopfl (iu21, gADFFiles%main)
      call kfopsc (iu21, 'POLTDDFT')
   endif

   if (master()) then
      idxlabel = 'dipole_x'
      call kfwrnr (iu21, idxlabel, gPoltddft%dipole(1:gPoltddft%dimecc,1), gPoltddft%dimecc, 1)
      idxlabel = 'dipole_y'
      call kfwrnr (iu21, idxlabel, gPoltddft%dipole(1:gPoltddft%dimecc,2), gPoltddft%dimecc, 1)
      idxlabel = 'dipole_z'
      call kfwrnr (iu21, idxlabel, gPoltddft%dipole(1:gPoltddft%dimecc,3), gPoltddft%dimecc, 1)
      idxlabel = 'magnetic_x'
      call kfwrnr (iu21, idxlabel, gPoltddft%magnetic(1:gPoltddft%dimecc,1), gPoltddft%dimecc, 1)
      idxlabel = 'magnetic_y'
      call kfwrnr (iu21, idxlabel, gPoltddft%magnetic(1:gPoltddft%dimecc,2), gPoltddft%dimecc, 1)
      idxlabel = 'magnetic_z'
      call kfwrnr (iu21, idxlabel, gPoltddft%magnetic(1:gPoltddft%dimecc,3), gPoltddft%dimecc, 1)
   endif

   allocate(intfitc(gDims%nsfos))
   intfitc = intfit

   allocate(gPoltddft%ranalysis(gPoltddft%dimecc))
   gPoltddft%ranalysis = 0.0_KREAL
   allocate(gPoltddft%ianalysis(gPoltddft%dimecc))
   gPoltddft%ianalysis = 0.0_KREAL
   allocate(gPoltddft%alpr(gPoltddft%npoints,3))
   gPoltddft%alpr = 0.0_KREAL
   allocate(gPoltddft%alpi(gPoltddft%npoints,3))
   gPoltddft%alpi = 0.0_KREAL
   allocate(gPoltddft%betr(gPoltddft%npoints,3))
   gPoltddft%betr = 0.0_KREAL
   allocate(gPoltddft%beti(gPoltddft%npoints,3))
   gPoltddft%beti = 0.0_KREAL

!start do npoints
   
   do ien = 1, gPoltddft%npoints

      call NewDistributedMatrix(Qr,REAL_MATRIX_TYPE,gDims%nsfos,gDims%nsfos,context,block)
      call NewDistributedMatrix(Qi,REAL_MATRIX_TYPE,gDims%nsfos,gDims%nsfos,context,block)     
      call NewDistributedMatrix(Dk,REAL_MATRIX_TYPE,gDims%nsfos,gDims%nsfos,context,block)
      call SetConstant(Qr, 0.0_KREAL)
      call SetConstant(Qi, 0.0_KREAL)
      call SetConstant(Dk, 0.0_KREAL)

      allocate(dc(gDims%nsfos,3))
      dc = zero_c

      allocate(dm(gDims%nsfos,3))
      dm = zero_c

      w = CMPLX(gPoltddft%wr(ien),gPoltddft%wi,kind=KREAL)

      ifi = gPoltddft%nkgrid 
      do i = 1, gPoltddft%nkgrid
          if (gPoltddft%intervals(i) >= gPoltddft%wr(ien)+gPoltddft%cutoff) then
             ifi = i
             exit
          endif
      enddo

      n = 0

      do i = 1, ifi
         
         if (gPoltddft%kdim(i) == 0) cycle

         tempc = (gPoltddft%intervals(i) + gPoltddft%intervals(i+1))/2
         sk = 2*(one_c/(tempc-w) + one_c/(tempc+w))
         skr = sk
         ski = AIMAG(sk)

!--------------------------------------------------
! analysis 2

         term1 = one_c/(tempc-w)
         term2 = one_c/(tempc+w)

         tk = 2*(term1-term2)
         tkr = tk
         tki = AIMAG(tk)

         tk = CMPLX(-tkr,tki)

!--------------------------------------------------

         allocate(sharedAk(gDims%nsfos,gPoltddft%kdim(i)))
         sharedAk = 0.0_KREAL
         call SetConstant(Dk, 0.0_KREAL)

         do k = n+1, n+gPoltddft%kdim(i)
           kstr = ' '
           call csputi (kstr, k)
           call kfread(iu63, 'Ak'//trim(kstr), sharedAk(1:gDims%nsfos,k-n))
         end do

         n = n + gPoltddft%kdim(i)

         call CalcFromFullMatrixProduct(Dk, sharedAk, sharedAk, .false., .true.)

         deallocate(sharedAk)

         call MultAdd(Qr, Dk, skr)
         call MultAdd(Qi, Dk, ski)

         dc(:,1) = dc(:,1) + sk*(gPoltddft%Vx(:,i))
         dc(:,2) = dc(:,2) + sk*(gPoltddft%Vy(:,i))
         dc(:,3) = dc(:,3) + sk*(gPoltddft%Vz(:,i))

!--------------------------------------------------
! analysis 2

         dm(:,1) = dm(:,1) + tk*(gPoltddft%Vmx(:,i))
         dm(:,2) = dm(:,2) + tk*(gPoltddft%Vmy(:,i))
         dm(:,3) = dm(:,3) + tk*(gPoltddft%Vmz(:,i))

!--------------------------------------------------

      enddo

      call DeleteDistributedMatrix(Dk)

      call NewDistributedMatrix(L,REAL_MATRIX_TYPE,gDims%nsfos,gDims%nsfos,context,block)
      call SetConstant(L, 0.0_KREAL)
      call LoadDistributedMatrix (L, iu61, 1, 'L')

      call CalcFromMatrixProduct(Mr, Qr, L)
      call CalcFromMatrixProduct(Mi, Qi, L)

      call DeleteDistributedMatrix(Qr)
      call DeleteDistributedMatrix(Qi)

      call NewDistributedMatrix(S,REAL_MATRIX_TYPE,gDims%nsfos,gDims%nsfos,context,block)
      call SetConstant(S, 0.0_KREAL)
      call LoadDistributedMatrix (S, iu61, 1, 'S')
      call InplaceAdd(Mr, S)
      call DeleteDistributedMatrix(S)

      call NewDistributedMatrix(C, Mr, Mi)
      
      call DeleteDistributedMatrix(Mr)
      call DeleteDistributedMatrix(Mi)

      do m = 1, 3    

         allocate(dp(gDims%nsfos,3))
         dp = zero_c

         dp(:,1) = dc(:,m)
         dp(:,2) = intfitc(:)

!-----------------------------
! analysis 2
         dp(:,3) = dm(:,m)
!-----------------------------

         info = 0 
         call SolveLinearSystemComplex(C, vectors=dp, info=info)
         if (info /= 0) then
            write(iuout,*)'SolveLinearSystemComplex returned info ', info
            call Print(C, iuout, 'The matrix in question')
            call stopit('SolveLinearSystemComplex returned info /= 0')
         end if

         vala = sum(intfitc*dp(:,1))
         valb = sum(intfitc*dp(:,2))

!-----------------------------------
! analysis 2
         vald = sum(intfitc*dp(:,3))
!-----------------------------------

         valc = vala/valb

         dc(:,m) = dp(:,1) - valc*dp(:,2)

         allocate(br(gDims%nsfos))
         br = 0.0_KREAL
         allocate(bi(gDims%nsfos))
         bi = 0.0_KREAL

         br = dc(:,m)
         bi = AIMAG(dc(:,m))

!----------------------------------
! analysis 2
        
         valc = vald/valb
         dm(:,m) = dp(:,3) - valc*dp(:,2)
         deallocate(dp)

         allocate(gr(gDims%nsfos))
         gr = 0.0_KREAL
         allocate(gi(gDims%nsfos))
         gi = 0.0_KREAL

         gr = dm(:,m)
         gi = AIMAG(dm(:,m))



!----------

!Save on TAPE bi to plot the induced density

         if (master()) then
            if (m == 1) then
               idxlabel = 'bx' 
            else if ( m == 2) then 
               idxlabel = 'by'
            else if ( m == 3) then
               idxlabel = 'bz'
            endif
            call csaddi ( idxlabel, '_', ien)
            call kfwrnr (iu21, idxlabel, bi, gDims%nsfos, 1)         
         endif

! save analysis 2
    
         if (master()) then
            if (m == 1) then
               idxlabel = 'qx'
            else if ( m == 2) then
               idxlabel = 'qy'
            else if ( m == 3) then
               idxlabel = 'qz'
            endif
            call csaddi ( idxlabel, '_', ien)
            call kfwrnr (iu21, idxlabel, gi, gDims%nsfos, 1)
         endif

         deallocate(gr)
         deallocate(gi)

!---------

         allocate(vectr(gDims%nsfos))
         allocate(vecti(gDims%nsfos))
         vectr = 0.0_KREAL
         vecti = 0.0_KREAL

         call MatrixVectorProduct(L, br, vectr)
         call MatrixVectorProduct(L, bi, vecti)
         deallocate(br)
         deallocate(bi)

         gPoltddft%ianalysis = 0.0_KREAL

         n = 0
         ki = 1
         do i = 1, ifi

            if (gPoltddft%kdim(i) == 0) cycle
            kd = gPoltddft%kdim(i)

            tempc = (gPoltddft%intervals(i) + gPoltddft%intervals(i+1))/2
            term1 = one_c/(tempc-w)
            term2 = one_c/(tempc+w)

            sk = 2*(term1+term2)
            skr = sk
            ski = AIMAG(sk)

            allocate(vecAr(kd))
            allocate(vecAi(kd))
            allocate(testr(kd))
            allocate(testi(kd))
            allocate(sharedAk(gDims%nsfos,kd))

            vecAr = 0.0_KREAL
            vecAi = 0.0_KREAL
            testr = 0.0_KREAL
            testi = 0.0_KREAL
            sharedAk = 0.0_KREAL

            do k = n+1, n+kd
              kstr = ' '
              call csputi (kstr, k)
              call kfread(iu63, 'Ak'//trim(kstr), sharedAk(1:gDims%nsfos,k-n))
            end do
            n = n + gPoltddft%kdim(i) 

            vecAr = matmul(vectr,sharedAk)
            vecAi = matmul(vecti,sharedAk)

            deallocate(sharedAk)

!-------------------

            kiEnd = ki+kd-1
            vecAr = -vecAr + gPoltddft%dipole(ki:kiEnd,m)
            testr =  skr * vecAr + ski * vecAi
            testi = -skr * vecAi + ski * vecAr
   
            gPoltddft%alpr(ien,m) = gPoltddft%alpr(ien,m) + sum(gPoltddft%dipole(ki:kiEnd,m)*testr)
            gPoltddft%alpi(ien,m) = gPoltddft%alpi(ien,m) + sum(gPoltddft%dipole(ki:kiEnd,m)*testi)

!------------------
!For analysis

!           gPoltddft%ranalysis(ki:kiEnd) = testr ! PIER
            gPoltddft%ianalysis(ki:kiEnd) = testi

!-------------------
!CD Spectra

            tk = 2*(term1-term2)
            tkr = tk
            tki = AIMAG(tk)

            testr =  tkr * vecAr + tki * vecAi
            testi = -tkr * vecAi + tki * vecAr
            gPoltddft%ranalysis(ki:kiEnd) = testi ! PIER

            gPoltddft%betr(ien,m) = gPoltddft%betr(ien,m) + sum(gPoltddft%magnetic(ki:kiEnd,m)*testr)
            gPoltddft%beti(ien,m) = gPoltddft%beti(ien,m) + sum(gPoltddft%magnetic(ki:kiEnd,m)*testi)            

!---------------
            deallocate(vecAr)
            deallocate(vecAi)
            deallocate(testr)
            deallocate(testi)
            ki = ki + kd
 
         enddo
         deallocate(vectr)
         deallocate(vecti)

         if (master()) then
            if (m == 1) then
               idxlabel = 'TCMx'
            else if ( m == 2) then
               idxlabel = 'TCMy'
            else if ( m == 3) then
               idxlabel = 'TCMz'
            endif
            call csaddi (idxlabel, '_', ien)
            call kfwrnr (iu21, idxlabel, gPoltddft%ianalysis, gPoltddft%dimecc, 1)
         endif

         if (master()) then
            if (m == 1) then
               idxlabel = 'TCMxR'
            else if ( m == 2) then
               idxlabel = 'TCMyR'
            else if ( m == 3) then
               idxlabel = 'TCMzR'
            endif
            call csaddi (idxlabel, '_', ien)
            call kfwrnr (iu21, idxlabel, gPoltddft%ranalysis, gPoltddft%dimecc, 1)
         endif

      enddo

      call DeleteDistributedMatrix(L)
      call DeleteDistributedMatrix(C)

      deallocate(dc)
      
!-------------------
! analysis 2

      deallocate(dm)

!-------------------

   enddo

   call ppNodebarrier

   if (master()) then
      call kfclfl(iu21)      

      call PrintandWriteKFPolTDDFT
   endif

   call kfclfl(iu61)
   call kfclfl(iu63)  
    
   deallocate(gPoltddft%ranalysis)
   deallocate(gPoltddft%ianalysis)   
   
   deallocate(intfit)
   deallocate(intfitc)

   call adf_end_blacs(context)

   call timere ('pp')
   call timere ('Alp')

End subroutine CalcAlp
