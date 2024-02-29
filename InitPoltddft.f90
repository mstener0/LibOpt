subroutine InitPoltddft(fragments)

   use SymUser
   use ResponseInterfaces
   use ExcitationInterfaces
   use KF
   use DimensionsModule
   use GeometryTypeModule
   use FragmentsArrayTypeModule
   use NucleiModule
   use FitsymModule
   use SymptrcModule
   use OptionModule
   use ExcitationModule
   use SpinorbitModule
   use ExcitOpenMod
   use CharsaModule
   use SpinsuffixModule
   use PoltddftModule
   use ADFGlobalInputModule
   use ADFFilesModule

   use ResponseCommon

   use RestartdataModule
   use ppCommonsModule

   use Vartypes
   implicit none

   type(FragmentsArrayType), intent(inout) :: fragments

   real(KREAL), parameter :: TWO   = 2.0_KREAL
   real(KREAL), parameter :: small = 1e-6_KREAL

   integer(KINT) :: noccmx(2), nunmax(2), isym, irsym
   integer(KINT) :: ispin, isos, nsos

   integer(KINT) :: iexist

   real(KREAL)   :: halffull, froc(2*gDims%naos)

   integer(KINT) :: iu, isetat, npeqbk, nsatbk, natsbk(gDims%nsetat), noatbk(gDims%nnuc),    &
          notypsbk(gDims%nsetat), nspin0

   logical                ::  lks, tda, nonaufbau, LEND63t


   LEND63t = .false.
   iexist = 0
   if (kfexfl('TAPE63')) iexist = 1

   call ppNodeBarrier
   call ppcbi(iexist,'iexst')

   if (iexist /= 0) then
      if ( .not. gRestartdata%restrt) call msg('WARNING: you should restart from TAPE21 generated from the same &
                                               simulation which generated TAPE63 in order to guarantee orbital phase &
                                               consistency') !there were a stopit here
      call kfopflpp(iu, 'TAPE63')
      call kfclfl(iu)
   endif  

   if (gFitsym%lnosft) then
      npeqbk = gDims%npeq
      nsatbk = gDims%nsetat
      do isetat = 1, gDims%nsetat
         natsbk(isetat) = gNuclei%nratst(isetat)
      end do
      noatbk(1:gDims%nnuc)     = gNuclei%noat(1:gDims%nnuc)
      notypsbk(1:gDims%nsetat) = gNuclei%notyps(1:gDims%nsetat)
      call kfopfl (iu, gADFFiles%main)
      call kfopsc (iu, 'SymFit')
      call kfread (iu, 'npeq', gDims%npeq)
      call kfread (iu, 'nsetat', gDims%nsetat)
      call kfrdni (iu, 'nratst', gNuclei%nratst, gDims%nsetat, 1)
      call kfrdni (iu, 'noat', gNuclei%noat, gDims%nnuc, 1)
      call kfrdni (iu, 'notyps', gNuclei%notyps, gDims%nsetat, 1)
      call kfclsc (iu)
      call kfclfl (iu)
   end if

!  --------------------------------------------------------
!  Determination of number of occupied and virtual orbitals
!  (noccmx and nunmax)
!  Hack!  in case of KSSPECTRUM:
!     the global arrays norboc or nroroc are modified such 
!     that they also work for partially occupied orbitals
!  --------------------------------------------------------

   noccmx = 0
   nunmax = 0
   call gInput%Get ('KSSpectrum', lks)
   call gInput%Get ('TDA', tda)
   tda = tda .or. lks
   nonaufbau = lks .or. gInput%IsValid ('ALLOW NONAUFBAU')

   if (gModel%ioprel<10) then

      if (gModel%nspin==1) gSymptrc%norboc(2,:) = 0

      if (nonaufbau) then
         halffull = (two - gOption%iun)/two
         call kfopfl (iu, gADFFiles%main)
         do ispin = 1, gModel%nspin
            do isym = 1, gDims%nsym
               gSymptrc%norboc(ispin,isym) = 0
               nsos = gSymptrc%norb(isym)
               call kfopsc (iu, gCharsa%bb(isym))
               call kfrdnr (iu, 'froc'//gSpinsuffix%sspin(ispin), froc, nsos, 1)
               do isos = 1, nsos
                  if (froc(isos)+small>halffull) then
                     gSymptrc%norboc(ispin,isym) = gSymptrc%norboc(ispin,isym) + 1
                  end if
               end do
            end do
         end do
         call kfclfl (iu)
      end if
      do ispin = 1, gModel%nspin
         do isym = 1, gDims%nsym
           noccmx(ispin) = noccmx(ispin) + gSymptrc%norboc(ispin,isym)
         end do
      end do
      if (gModel%nspin==1) noccmx(2) = 0

      if (abs(gModel%qelec-nint(gModel%qelec))>small .and. .not. nonaufbau) then
         call stopit ('A fractional number of electrons detected in the subroutine Symaexcit')
      end if
      if ((noccmx(1)*(2-gOption%iun)+noccmx(2)*gOption%iun-nint(gModel%qelec))/=0) then
         write (iuout,9000) noccmx(1), noccmx(2), gOption%iun, gModel%qelec
         if (.not.nonaufbau) call stopit                                                          &
         ('Incorrect number of occupied orbitals calculated. Probably no aufbau solution, which is needed.')
      end if

      if (gModel%nspin==2) then
         nunmax = gDims%nsot - noccmx
      else
         nunmax(1) = gDims%nsot - noccmx(1)
      endif

   else
      if (nonaufbau) then
         call kfopfl (iu, gADFFiles%main)
         do irsym = 1, gSpinorbit%nrsym
            gSpinorbit%nroroc(1,irsym) = 0
            nsos = gSpinorbit%nrorb(irsym)
            call kfopsc (iu, gSpinorbit%bbr(irsym))
            if (nsos>0) call kfrdnr (iu, 'froc'//gSpinsuffix%sspin(1), froc, nsos, 1)
            halffull = gSpinorbit%nrdim(irsym)/two
            do isos = 1, nsos
               if (froc(isos)+small>halffull) then
                  gSpinorbit%nroroc(1,irsym) = gSpinorbit%nroroc(1,irsym) + 1
               end if
            end do
         end do
         call kfclfl (iu)
      end if
      do irsym = 1, gSpinorbit%nrsym
         noccmx(1) = noccmx(1) + gSpinorbit%nroroc(1,irsym)*gSpinorbit%nrdim(irsym)
         nunmax(1) = nunmax(1) + (gSpinorbit%nrorb(irsym)-gSpinorbit%nroroc(1,irsym))              &
               *gSpinorbit%nrdim(irsym)
      end do
      noccmx = 2*noccmx
   end if

   if (gModel%lSpinFlipTDDFT) then
      nspin0 = nunmax(2)
      nunmax(2) = nunmax(1)
      nunmax(1) = nspin0
   end if

   call CalcIntegrals(noccmx, nunmax, iexist, fragments, LEND63t)
   if (.not.LEND63t) then
      call CalcVxyz()
      call msg ('MATRIXAK')
      call CalcAlp()

      deallocate (gPoltddft%skipnb) 
      deallocate (gPoltddft%noccmin, gPoltddft%nvirtmax, gPoltddft%kidx, gPoltddft%kdim)
      deallocate (gPoltddft%energy, gPoltddft%intervals, gPoltddft%wr)
      call DeleteSharedArray(gPoltddft%Vx, 'Vx')
      call DeleteSharedArray(gPoltddft%Vy, 'Vy')
      call DeleteSharedArray(gPoltddft%Vz, 'Vz')
      call DeleteSharedArray(gPoltddft%Vmx, 'Vmx')
      call DeleteSharedArray(gPoltddft%Vmy, 'Vmy')
      call DeleteSharedArray(gPoltddft%Vmz, 'Vmz')
      call DeleteSharedArray(gPoltddft%dipole,'dipole')
      call DeleteSharedArray(gPoltddft%dipvel,'dipvel')
      call DeleteSharedArray(gPoltddft%magnetic,'magnetic')
      deallocate(gPoltddft%alpr)
      deallocate(gPoltddft%alpi)
      deallocate(gPoltddft%betr)
      deallocate(gPoltddft%beti)
   end if

   call kfdlfl('TAPE58')
   !call kfdlfl('TAPE59')
   call kfdlfl('TAPE60')
   call kfdlfl('TAPE61')

   if (gFitsym%lnosft) then
      gDims%npeq   = npeqbk
      gDims%nsetat = nsatbk
      gNuclei%nratst(1:gDims%nsetat) = natsbk(1:gDims%nsetat)
      gNuclei%notyps(1:gDims%nsetat) = notypsbk(1:gDims%nsetat)
      gNuclei%noat(1:gDims%nnuc)     = noatbk(1:gDims%nnuc)
   end if

 9000 format(' noccmx, nocmx1, iun and qelec in Symaexcit are equal to',3I4,F15.10)
end subroutine InitPoltddft
