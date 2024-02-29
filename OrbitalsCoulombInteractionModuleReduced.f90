module OrbitalsCoulombInteractionModuleReduced


   use VarTypes
   use FragmentTypeModule   ! here just use only FragmentTypeModule and not even FragmentsArrayTypeModule
                            ! FragmentsArrayTypeModule includes FragmentTypeModule

   use AdfDensityEvaluatorModule
   use KF
   use ZlmFitModule
   use ZlmFitConfigModule
   use NewZlmFitConfigAdfModule
   use AssertionsModule
   use DimensionsModule
   use InputReaderModule
   use ftlStringModule
   use SharedArraysModule
   use PMatForSingleOrbitalModule
   use M_grid

   implicit none

   public :: CalcOrbitalsCoulombInteractionReduced
   public :: gOrbitalsCoulombInteractionReduced

   type arrayPotentialI
       real(KREAL), allocatable :: potentialI(:)
   end type arrayPotentialI

   type OrbitalsCoulombInteractionReducedType
      integer(KINT) :: icounter
      real(KREAL)   :: totalRhoI
   end type
   type(OrbitalsCoulombInteractionReducedType), save :: gOrbitalsCoulombInteractionReduced

   type(arrayPotentialI), allocatable, save :: giBlock(:)


contains

subroutine CalcOrbitalsCoulombInteractionReduced(fragment, iOrbital, jOrbital, correctionenergy) ! automatic acquisition

   implicit none

   type(fragmentType), intent(inout) :: fragment    ! fragment = fragments%activeFragment
   integer(KINT), intent(inout) :: iOrbital, jOrbital      
   real(KREAL), intent(inout) :: correctionenergy   ! for pass correction

   ! ============================================================================
   ! Calculate the electrostatic (coulomb) interaction energy between the density
   ! of two orbitals.
   ! ============================================================================

   type(AdfDensityEvaluatorType) :: orbitalEvaluatorI, orbitalEvaluatorJ
   integer :: iBlock
   integer :: iSpin = 1   ! can not do unrestricted, idem for fragment%nSpin, further on there is a stop      
   type(ZlmFitConfigType) :: zlmFitConfig
   type(ZlmFitType) :: zlmFit
   real(KREAL) :: interactionEnergy(fragment%nSpin), totalRhoJ(fragment%nSpin)
   real(KREAL), allocatable :: densityJ(:,:), densityI(:,:)
   real(KREAL), pointer :: pMatOrbitalI(:,:), pMatOrbitalJ(:,:)
   type(grid) :: G
   logical, parameter :: USE_EQUIVALENT_BLOCKS = .true.
   logical :: ievaluate
   real(KREAL), allocatable :: tmppotentialI(:)

   if (gOrbitalsCoulombInteractionReduced%icounter /= iOrbital) then  ! check counter for calculate potential and density of i-th occupied orbital one time
     ievaluate = .true.
     gOrbitalsCoulombInteractionReduced%icounter = iOrbital

     deallocate(giBlock) ! for calculate new potential and density of i-th occupied orbital one time, see previous routine 
   else
      ievaluate = .false.
   end if

   correctionenergy = 0.0_KREAL

   call TCAssert(fragment%nSpin==1, 'OrbitalsCoulombInteraction: Cannot do spin-unrestricted calculations (yet).')  ! necessary? already present?

   if (ievaluate) then  ! for calculate potential and density of new occupied orbital
     call NewSharedArray(pMatOrbitalI, 'pMatOrbitalI', [gDims%naosx,fragment%nSpin])
   end if

   call NewSharedArray(pMatOrbitalJ, 'pMatOrbitalJ', [gDims%naosx,fragment%nSpin])
 
   call create (G, USE_EQUIVALENT_BLOCKS)
   call NewZlmFitConfigAdf (zlmFitConfig, fragment, printZlmFitConfigInfo=.false.)
   call zlmFit%New(zlmFitConfig, fragment%geom%Coordinates(:,:,1), nint(fragment%geom%qtch,KINT), fragment%nSpin)
 
!   write(iuout,*) "***********************************"
!   write(iuout,*) ""
!   write(iuout,*) "Computing interaction energy between the density of orbital #"//ftlString(iOrbital)//' and orbital #'//ftlString(jOrbital)

   ! Compute the integral between the potential of orbital I and and the density of orbital J
   ! ========================================================================================

   if (ievaluate) then    ! for calculate potential and density of new occupied orbital

     call GetPMatForSingleOrbital(pMatOrbitalI, fragment, iOrbital)
     call orbitalEvaluatorI%New(fragment, pmatBas=pMatOrbitalI, onlyValence=.true.)
     call zlmFit%CalcFit (orbitalEvaluatorI)

     call GetPMatForSingleOrbital(pMatOrbitalJ, fragment, jOrbital)
     call orbitalEvaluatorJ%New(fragment, pmatBas=pMatOrbitalJ, onlyValence=.true.)

     gOrbitalsCoulombInteractionReduced%totalRhoI = 0.0_KREAL
     totalRhoJ = 0.0_KREAL
     interactionEnergy = 0.0_KREAL

   else   ! potential and density of occupied orbital is just calculated

     call GetPMatForSingleOrbital(pMatOrbitalJ, fragment, jOrbital)
     call orbitalEvaluatorJ%New(fragment, pmatBas=pMatOrbitalJ, onlyValence=.true.)

     totalRhoJ = 0.0_KREAL
     interactionEnergy = 0.0_KREAL

   end if

   call timers ('pp')

   if (ievaluate) then ! calculate and save potential and density of new occupied orbital

     allocate(giBlock(G%nblocks))
     allocate(tmppotentialI(G%npoints))
     allocate(densityJ(G%npoints,fragment%nSpin), densityI(G%npoints,fragment%nSpin))

     iblock_: do iBlock = 1, G%nblocks

       if (skip_block(G, iBlock)) cycle iblock_

       call get_block(G)

       ! Compute density of orbital J
       call orbitalEvaluatorJ%Calc(G%coord(:,:), densityJ)

       ! Compute density and potential of orbital I
       call zlmfit%CalcDenPot(G%coord(:,:), density=densityI, potential=tmppotentialI)
       giBlock(iBlock)%potentialI = tmppotentialI
       gOrbitalsCoulombInteractionReduced%totalRhoI =  gOrbitalsCoulombInteractionReduced%totalRhoI + sum(densityI(:,iSpin)*g%w)

       totalRhoJ(iSpin) = totalRhoJ(iSpin) + sum(densityJ(:,iSpin)*g%w)
       interactionEnergy = interactionEnergy + sum(giBlock(iBlock)%potentialI*densityJ(:,iSpin)*g%w)

     end do iblock_
     deallocate(tmppotentialI)
     deallocate(densityJ, densityI)

     call ppcbnr (gOrbitalsCoulombInteractionReduced%totalRhoI, iSpin, 'totalRhoI')  ! only restricted 

   else   ! potential and density of occupied orbital is just calculated

     allocate(densityJ(G%npoints,fragment%nSpin))

     iblock2_: do iBlock = 1, G%nblocks

       if (skip_block(G, iBlock)) cycle iblock2_

       call get_block(G)

       ! Compute density of orbital J
       call orbitalEvaluatorJ%Calc(G%coord(:,:), densityJ)

       totalRhoJ(iSpin) = totalRhoJ(iSpin) + sum(densityJ(:,iSpin)*g%w)
       interactionEnergy = interactionEnergy + sum(giBlock(iBlock)%potentialI*densityJ(:,iSpin)*g%w)

     end do iblock2_
     deallocate(densityJ)

   end if

   call rewind_grid (G)

   call ppcbnr (totalRhoJ, size(totalRhoJ), 'totalRhoJ')

   call ppcbnr (interactionEnergy, 1, 'interactionEnergy')
   call timere ('pp')    

!   write(iuout,*) ""
!   write(iuout,*) "Coulomb interaction energy [a.u.]: ", interactionEnergy
!   write(iuout,*) ""
!   write(iuout,*) "Integrated density of orbital #"//ftlString(iOrbital)//" (should be approximately 1.0): ", gOrbitalsCoulombInteractionReduced%totalRhoI
!   write(iuout,*) "Integrated density of orbital #"//ftlString(jOrbital)//" (should be approximately 1.0): ", totalRhoJ
!   write(iuout,*) ""

   correctionenergy = interactionEnergy(1)   ! only restricted, for pass value

   if (ievaluate) then ! only for new occupied orbital
 
     call orbitalEvaluatorJ%Delete()
     call orbitalEvaluatorI%Delete()
     call DeleteSharedArray(pMatOrbitalI, 'pMatOrbitalI')
 
   else ! if occupied orbital is just calculated
 
     call orbitalEvaluatorJ%Delete()
 
   end if

   call DeleteSharedArray(pMatOrbitalJ, 'pMatOrbitalJ')

   call zlmFit%Delete()
   call zlmFitConfig%Delete ()
   call delete(G)

end subroutine CalcOrbitalsCoulombInteractionReduced

end module OrbitalsCoulombInteractionModuleReduced
