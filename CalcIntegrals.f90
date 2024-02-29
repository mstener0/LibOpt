subroutine CalcIntegrals(noccmx,nunmax, iexist, fragments, LEND63t)

   use KF
   use MasterModule
   use DimensionsModule
   use NucleiModule
   use DistributedMatrixModule
   use ScaMatrixArrayModule
   use SharedArraysModule
   use adf_blacs
   use ppCommonsModule
   use Vartypes

   use ResponseInterfaces,except=>IrredMatEltsMerge
   use ResponseModule
   use SymUser
   use SymptrcModule
   use SymmetrydataModule
   use ExcitOpenMod
   use FragmentsArrayTypeModule

   use ExcitModVars

   use ModelDataModule
   use ExcitSOVars
   use ExcitationModule
   use SymmetryInfo
   use SymrespstanModule

   implicit none 

   integer(KINT), intent(in)  :: noccmx(2), nunmax(2), iexist
   type(FragmentsArrayType), intent(inout) :: fragments

   logical, intent(inout) :: LEND63t

   logical, external :: master

   integer(KINT), allocatable :: nrocc(:,:), nrvirt(:,:)
   integer(KINT), allocatable :: nsymoc(:,:), nsymun(:,:)
   real(KREAL), allocatable   :: epsocc(:,:), epsun(:,:)
   integer(KINT), allocatable :: nspnoc(:,:), nspnun(:,:)
   logical, allocatable       :: selectocc(:,:), selectun(:,:)

   integer(KINT)   :: noccmax, nunomax, nds

   logical                    :: lmodify, lselect
   integer(KINT)              :: nocc2do(20), nvirt2do (20)
   integer(KINT), allocatable :: iocc2do(:, :), ivirt2do (:,:)

   logical                    :: lrep2do(20,MAXDIM,2)

   integer(KINT)   :: irrep
   integer(KINT)   :: i, idim
   integer(KINT)   :: nsymr, isym, i2, i3
   integer(KINT), allocatable :: irrepop(:)

!  copy from InitResponseSymmetry to calculate nrocc and nrvirt

   call timers ('Integrals')


   gExcitSOVars%nsymt = 0
   if (gModel%ioprel<10) then
      gExcitSOVars%nsymt = gDims%nsym
   else
      do isym = 1, gExcitSOVars%nrsym
         do idim = 1, gExcitSOVars%nrdim(isym)
            gExcitSOVars%nsymt = gExcitSOVars%nsymt + 1
         end do
      end do
      i     = 0
      gExcitSOVars%ndimt = 0
      do isym = 1, gExcitSOVars%nrsym
         do idim = 1, gExcitSOVars%nrdim(isym)
            i = i + 1
            if (idim==1) gExcitSOVars%ndimt(i) = gExcitSOVars%nrdim(isym)
         end do
      end do
   end if

   if (gModel%ioprel<10) nsymr = gDims%nsym
   if (gModel%ioprel>=10) nsymr = gExcitSOVars%nsymt

   noccmax = max0(noccmx(1),noccmx(2))
   nunomax = max0(nunmax(1),nunmax(2))

   nds = 1
   if (gModel%nspin==2) nds = 2
   gExcitOpenMod%nds = nds
   gExcitOpenMod%noccmax = noccmax
   gExcitOpenMod%nunomax = nunomax

   allocate (nrocc(nsymr,nds))
   allocate (nrvirt(nsymr,nds))
   allocate (nsymoc(noccmax,nds))
   allocate (nsymun(nunomax,nds))

   nrocc  = 0
   nrvirt = 0 
   nsymoc = 0
   nsymun = 0

   allocate (epsocc(noccmax,nds), epsun(nunomax,nds))
   allocate (nspnoc(noccmax,nds), nspnun(nunomax,nds))
   allocate (selectocc(noccmax,nds), selectun(nunomax,nds))
   allocate (gExcit%iorboc(noccmax,nds), gExcit%norbun(nunomax,nds))
   allocate (iocc2do(gDims%naos, 20), ivirt2do (gDims%naos, 20))

   epsocc = 0.0_KREAL
   epsun = 0.0_KREAL

   call SetModifyExcitation (lmodify, lselect, nocc2do, nvirt2do, iocc2do, ivirt2do )

   call ReadEpsilons (noccmax, nunomax, epsocc, epsun, nsymoc, nsymun, &
                      nspnoc, nspnun, gExcit%iorboc,     &
                      gExcit%norbun, lmodify, lselect, nocc2do, nvirt2do, iocc2do, ivirt2do, &
                      selectocc, selectun )

   deallocate (gExcit%iorboc, gExcit%norbun)
   deallocate (nspnoc, nspnun)
   deallocate (iocc2do, ivirt2do )

   do isym = 1, nsymr
      do i2 = 1, noccmax
         if (nsymoc(i2,1)==isym) nrocc(isym, 1) = nrocc(isym, 1) + 1
      enddo
      do i3 = 1, nunomax
         if (nsymun(i3,1)==isym) nrvirt(isym, 1) = nrvirt(isym, 1) + 1
      end do
   end do

   call InitResponseSymmetry(nsymoc, nsymun, noccmx, nunmax, nrocc, nrvirt)

   call SetPoltddftValues(lrep2do)
   allocate(irrepop(gSymmetryInfo%nrep))
   irrepop = 0

   do irrep = 1, gSymmetryInfo%nrep
      if(lrep2do(irrep,1,1)) irrepop(irrep) = irrep
   enddo

   deallocate (nsymoc,nsymun)

   call CalcMatrixL ()

   call CalcPairfit()

   call ppNodeBarrier

   call CreateKGrid(noccmx, nunmax, nsymr, nrocc, nrvirt, epsocc, epsun, irrepop, fragments)
   call ppNodeBarrier
write(iuout,*) "In CalcINtegrals after CreateKGrid PIER"
   if (iexist ==0) then
write(iuout,*) "In CalcINtegrals inside if for CalcAk PIER"
      call CalcAk(nsymr, nrocc, nrvirt, irrepop, LEND63t)
   endif
write(iuout,*) "In CalcINtegrals after if for CalcAk PIER"

   call ppNodeBarrier
write(iuout,*) "PIER 1"
   if (.not.LEND63t) then
write(iuout,*) "PIER 2"
      call CalcDip(nsymr, nrocc, nrvirt, irrepop, iexist)
write(iuout,*) "PIER 3"
      call ReorderKGrid(iexist)
write(iuout,*) "PIER 4"
   end if
   deallocate(irrepop)
   deallocate (nrocc, nrvirt)
   deallocate (epsocc, epsun)
   deallocate (selectocc, selectun)
   call timere ('Integrals')

write(iuout,*) "end of CalcINtegrals PIER"

end subroutine CalcIntegrals
