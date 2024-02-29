subroutine CalcMatrixL ()
!
   use KF
   use MasterModule
   use DimensionsModule
   use SCFUtilsModule
   use DistributedMatrixModule
   use adf_blacs
   use PoltddftModule
   use ppCommonsModule
!
   use Vartypes
   use ADFGlobalInputModule
   use ModelDataModule, only: gModel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
   use Vartypes
   use HFParallelizationModule
   use ModelDataModule
   use XCFunctionalDescriptorModule
   use RangeSeparatedCalculation
   use XCRangeSeparatedCalcDescriptor
   use fitfitint2RangeSeparatedModule
   use ADFFilesModule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   implicit none

   type(DistributedMatrixType) :: S, F, Z, eigvec, Ztransp, S2, F2, Frs, S3, Z2, Zrs

   real(KREAL)   :: eigval(gDims%nsfos)
   integer(KINT) :: context1, block
   integer(KINT) :: info
   integer(KINT) :: iu, iu65, iu69
   logical       :: HDA_fitted
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
   integer(KINT) :: iurs, iu777
   logical       :: IsRangeSeparated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer(KINT) :: iexist64

!  ---------------------------------
!  opem tape with fit integrals sfit
!  ---------------------------------

   call GetBLACSContextAndBlocksize(gDims%nsfos, context1, block)

   call timers ('pp')

   call kfopfl (iu, 'TAPE59')
   call ReadPacLowAOMatrixToDist(iu, 'FitFit%fitfit', F, gDims%nsfos, context1, block)
   call kfclfl (iu)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
   IsRangeSeparated = IsHFRangeSeparatedFunctional(gModel%xcDescriptor)

   if (IsRangeSeparated) then
      call kfopfl (iurs, 'TAPE77')
      call ReadPacLowAOMatrixToDist(iu, 'FitFitRS%fitfitRS', Frs, gDims%nsfos, context1, block)
      call kfclfl (iurs)
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(iuout,*) "matrix F elements PIER", F%rValues
!write(iuout,*) "matrix Frs elements PIER", Frs%rValues

   call kfopfl (iu, 'TAPE60')
   call ReadPacUprAOMatrixToDist(iu, 'MatrixZ%Z', Z, gDims%nsfos, context1, block)
   call kfclfl (iu)
!!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
if (IsRangeSeparated) call NewDistributedMatrix(Z2, Z, copyData=.true.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call InplaceAdd(Z, F)
   call Scale(Z,gPoltddft%lambda)

!!!!!!!!!!!!!!!!PIERPAOLO
if (IsRangeSeparated) then
   call kfopfl (iu777, 'TAPE777')
   call ReadPacUprAOMatrixToDist(iu777, 'MatrixZrs%Zrs', Zrs, gDims%nsfos, context1, block)
   call kfclfl (iu777)
   call InplaceAdd(Z2, Frs)
   call InplaceSubtract(Z2 , Zrs )
   call Scale(Z2,gPoltddft%lambda)
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call DeleteDistributedMatrix(F)

   call kfopfl (iu, 'TAPE58')
   call ReadPacLowAOMatrixToDist(iu, 'Fitsfit%sfit', S, gDims%nsfos, context1, block)
   call kfclfl (iu)

   call kfcrfl (iu, 'TAPE61')
   call kfcrsc (iu, 'Matrix')

   call StoreDistributedMatrix (S, iu, 1, 'S')

   call gInput%Get ('PolTDDFT%HDA_fitted', HDA_fitted)

   if ((HDA_fitted).and.(.not. gModel%lhybrid)) call stopit('HDA_fitted support only hybrid functionals')

   if (kfexfl('TAPE64')) then
      iexist64 = 1
   end if
 
   call ppNodeBarrier
   call ppcbi(iexist64,'iexst64')
   if ((iexist64 ==1).and.(HDA_fitted)) then
      HDA_fitted = .false.
   end if

   if (HDA_fitted) then
!write(iuout,*) "creates S2 PIER"
      call NewDistributedMatrix(S2, S, copyData=.true.)
      call NewDistributedMatrix(S3, S, copyData=.true.)
   end if

   info = 0 
   if (context1 == CONTEXT_DUPLICATED_MATRIX) then
      call dposv('L',gDims%nsfos,gDims%nsfos,S%rValues,gDims%nsfos,Z%rValues,gDims%nsfos,info)
#ifdef D_SCALAPACK
   else if (context1 /= CONTEXT_PROC_OUTSIDE_CTXT) then
      call pdposv('L',gDims%nsfos,gDims%nsfos,S%rValues,1,1,S%desc,Z%rValues,1,1,Z%desc,info)
#endif
   endif
   if (info /= 0) then
      write(iuout,*) 'dposv info', info
      call LoadDistributedMatrix (S, iu, 1, 'S')
      call Diagonalize(S, eigval, eigvec, info)
      write(iuout,'(A/(10E14.5))')'Eigenvalues of the fit overlap matrix',eigval
      call stopit('CalcMatrixL: [p]dposv returned non-zero status')
   end if

   call StoreDistributedMatrix (Z, iu, 1, 'L') 

   if ((HDA_fitted).and.(.not.IsRangeSeparated)) then
      call NewDistributedMatrix(Ztransp, Z, copyData=.false.)
      call SetConstant(Ztransp, 0.0_KREAL)
      call InplaceAdd(Ztransp, Z, .true.)

!write(iuout,*) "inside .not.IsRangeSeparated PIER"

      if (context1 == CONTEXT_DUPLICATED_MATRIX) then
         call dposv('L',gDims%nsfos,gDims%nsfos,S2%rValues,gDims%nsfos,Ztransp%rValues,gDims%nsfos,info)
#ifdef D_SCALAPACK
      else if (context1 /= CONTEXT_PROC_OUTSIDE_CTXT) then
         call pdposv('L',gDims%nsfos,gDims%nsfos,S2%rValues,1,1,S2%desc,Ztransp%rValues,1,1,Ztransp%desc,info)
#endif
      endif

      call DeleteDistributedMatrix(S2)
      call kfcrfl(iu65,'TAPE65')
      call kfopfl(iu65,'TAPE65')
      call kfcrsc(iu65,'Q')
      call StoreDistributedMatrix (Ztransp, iu65, 1, 'Q')
      call kfclfl(iu65)

      call DeleteDistributedMatrix(Ztransp)

   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
   if (IsRangeSeparated) then

!write(iuout,*) "in IsRangeSeparated PIER"

   info = 0
   if (context1 == CONTEXT_DUPLICATED_MATRIX) then
      call dposv('L',gDims%nsfos,gDims%nsfos,S2%rValues,gDims%nsfos,Z2%rValues,gDims%nsfos,info)
#ifdef D_SCALAPACK
   else if (context1 /= CONTEXT_PROC_OUTSIDE_CTXT) then
      call pdposv('L',gDims%nsfos,gDims%nsfos,S2%rValues,1,1,S2%desc,Z2%rValues,1,1,Frs%desc,info)
#endif
   endif
   if (info /= 0) then
      write(iuout,*) 'dposv info', info
      call LoadDistributedMatrix (S, iu, 1, 'S')
      call Diagonalize(S, eigval, eigvec, info)
      write(iuout,'(A/(10E14.5))')'Eigenvalues of the fit overlap matrix',eigval
      call stopit('CalcMatrixL: [p]dposv returned non-zero status')
   end if

      call NewDistributedMatrix(Ztransp, Z2, copyData=.false.)
      call SetConstant(Ztransp, 0.0_KREAL)
      call InplaceAdd(Ztransp, Z2, .true.)

      if (context1 == CONTEXT_DUPLICATED_MATRIX) then
         call dposv('L',gDims%nsfos,gDims%nsfos,S3%rValues,gDims%nsfos,Ztransp%rValues,gDims%nsfos,info)
#ifdef D_SCALAPACK
      else if (context1 /= CONTEXT_PROC_OUTSIDE_CTXT) then
         call pdposv('L',gDims%nsfos,gDims%nsfos,S3%rValues,1,1,S3%desc,Ztransp%rValues,1,1,Ztransp%desc,info)
#endif
      endif

      call DeleteDistributedMatrix(S2)
      call DeleteDistributedMatrix(S3)
      call kfcrfl(iu65,'TAPE65')
      call kfopfl(iu65,'TAPE65')
      call kfcrsc(iu65,'Q')
      call StoreDistributedMatrix (Ztransp, iu65, 1, 'Q')
      call kfclfl(iu65)

      call DeleteDistributedMatrix(Ztransp)
       call DeleteDistributedMatrix(Frs)

   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call kfclfl (iu)

   call DeleteDistributedMatrix(S)
   call DeleteDistributedMatrix(Z)

   call timere ('pp')   

   call adf_end_blacs(context1)

!  ======================================================================

end subroutine CalcMatrixL

