subroutine SetPoltddftValues(lrep2dopoltddft)

!  ============================================================================
!  purpose: set defaults for parameters in POLTDDFT key block in ADF
!           INPUT file. The values are stored in a common block
!           poltddft.
!           Other variables are initialized. Information is read from
!           the input file. These variables are then given back to 
!           the calling program.
!
!  Input:   NONE
!
!  Output:
!           lrep2doexact    - what irrep, column will be calculated with
!                             the exact method (full diagonalization)
!           lcauchy         - Cauchy coeffiecients are calculated
!
!  called from: CalcIntegrals
!  ============================================================================

   use Vartypes
   use SymUser
   use ExcitationInterfaces,except=>SetValues
   use ModelDataModule
   use RespgeneralModule
   use ExcitationModule
   use CDSpectrum
   use MCDCommon
   use MCDInput
   use SymmetryInfo
   use PoltddftModule
   use ADFGlobalInputModule

   use Vartypes
   implicit none

   logical, external :: contns
   integer(KINT)     :: i
   real(KREAL)       :: step !backward_extention !PIERPAOLO

   logical, intent(out) :: lrep2dopoltddft(LNSYM,MAXDIM,2)

   real(KREAL) :: eVha
   real(KREAL), allocatable  :: rRange(:)
   real(KREAL), external :: convrs

   character(LCHARS), allocatable :: lines(:)

   call timers ('SetPoltddftValues')

   eVha = convrs('ENERGY','HARTREE,EV')

   lrep2dopoltddft   = .TRUE.

   call gInput%Get ('PolTDDFT%Velocity', gPoltddft%dv)
   call gInput%Get ('PolTDDFT%Lambda', gPoltddft%lambda)
   call gInput%Get ('PolTDDFT%Lifetime', gPoltddft%wi)
   call gInput%Get ('PolTDDFT%CutOff', gPoltddft%cutoff)
   call gInput%Get ('PolTDDFT%KGrid', gPoltddft%eVgrid)
   call gInput%Get ('PolTDDFT%NGrid', gPoltddft%nkgrid)
   call gInput%Get ('PolTDDFT%NFreq', gPoltddft%npoints)
   call gInput%Get ('PolTDDFT%FreqRange', rRange)
   if (size(rRange)/=2) call stopit('Expected lower and upper bound for FreqRange.')

   if (gInput%IsPresent ('POLTDDFT%Irrep')) then
      lrep2dopoltddft    = .FALSE.
      call gInput%GetFreeBlock ('POLTDDFT%IRREP', lines)
      call ReadRep2Do (.FALSE., lines, size(lines), lrep2dopoltddft)
   end if

   gPoltddft%cutoff = gPoltddft%cutoff/eVha
   gPoltddft%wi = gPoltddft%wi/eVha

   gPoltddft%eVgrid = gPoltddft%eVgrid/eVha
   step = gPoltddft%eVgrid/gPoltddft%nkgrid
!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
   !!backward_extention = 5.0d0/eVha
   !gPoltddft%nkgrid = gPoltddft%nkgrid + (backward_extention/step)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate (gPoltddft%intervals(gPoltddft%nkgrid+1))
   !gPoltddft%intervals(1) = -backward_extention !PIERPAOLO 0.0d0
   gPoltddft%intervals(1) = 0.0d0
   do i = 1, gPoltddft%nkgrid
      gPoltddft%intervals(i+1) = gPoltddft%intervals(i) + step
   enddo

   rRange = rRange/eVha
   allocate (gPoltddft%wr(gPoltddft%npoints+1))
   step = (rRange(2)-rRange(1))/gPoltddft%npoints
   gPoltddft%wr(1) = rRange(1)
   do i = 1, gPoltddft%npoints
      gPoltddft%wr(i+1) = gPoltddft%wr(i) + step
   enddo

   call timere ('SetPoltddftValues')

end subroutine SetPoltddftValues

