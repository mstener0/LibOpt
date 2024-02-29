Subroutine CalcDip(nsymr, nrocc, nrvirt, irrepop, iexist)

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
   use SharedArraysUtilModule
   use adf_blacs
   use SCFUtilsModule
   use ADFFilesModule

   use RestartdataModule

   use SymrespstanModule, only: gSymrespstan

   use PairsModule

   use MasterModule

   implicit none 

   integer(KINT), intent(in)  :: nsymr, irrepop(gSymmetryInfo%nrep)
   integer(KINT), intent(in)  :: nrocc(gExcitSOVars%nsymt,gExcitOpenMod%nds),&
                                 nrvirt(gExcitSOVars%nsymt,gExcitOpenMod%nds)

   integer(KINT) ,intent(in)  :: iexist

   integer(KINT)              :: i, minsize, irrep
   integer(KINT)              :: nocc, nvirt, nmo
   integer(KINT)              :: isym

   integer(KINT)              :: block, context1, iexist21 

   integer(KINT)              :: iu21, iu15

   integer(KINT)              :: nbas(nsymr), nsym
   character(len=LCHARS)      :: symrep(nsymr)

   integer(KINT), allocatable :: npart(:)

   real(KREAL), allocatable   :: tmp(:)
   real(KREAL), pointer       :: eigen(:,:)

   type(DistributedMatrixType):: Eocc(nsymr), Evirt(nsymr)

   call timers ('dipole')
   iexist21 = 0
   if (kfexfl('TAPE21')) iexist21 = 1

   if ((iexist.eq.1).and.(iexist21.eq.1)) then
      call kfopfl(iu21, gRestartdata%rsfile)
   end if
   if ((iexist.eq.1).and.(iexist21.eq.0)) then
      call kfopfl(iu21, gADFFiles%main)
   end if
   if (iexist.eq.0) then
      call kfopfl(iu21, gADFFiles%main)
   end if

   call kfread (iu21, 'Symmetry%nsym', nsym)      
   call kfread (iu21, 'symlab', symrep)
   call kfread (iu21, 'nfcn', nbas)
   minsize = min(minval(nrocc(1:nsymr, 1)),minval(nrvirt(1:nsymr, 1)))
   call GetBLACSContextAndBlocksize(minsize, context1, block)

   call timers('pp')
   do isym = 1, nsymr
      nocc = nrocc(isym, 1)
      nvirt = nrvirt(isym, 1)
      nmo = nocc+nvirt
      if (nmo == 0) cycle
      allocate(npart(nbas(isym)),tmp(nbas(isym)))
      npart = 0 
      tmp = 0.0_KREAL
      call kfopsc(iu21, symrep(isym))
      call kfread(iu21, 'npart', npart)
      call kfopvr(iu21, 'Eigen-Bas_A')
      ! Read eigenvectors and split them into occupied and virtuals, and expand to full naos

      call NewSharedArray(eigen,'eigen',(/gDims%naos,nmo/))
      if (AmIOwner(eigen)) then
         eigen = 0.0_KREAL
         do i = 1, nmo
            call kfread(iu21, '%', tmp)
            eigen(npart(1:nbas(isym)),i) = tmp(1:nbas(isym))
         end do
      end if
      call ppNodeBarrier
      if (nocc > 0) then
         call NewDistributedMatrix(Eocc(isym), eigen(:,1:nocc), context1, block)
      end if
      if (nvirt > 0) then
         call NewDistributedMatrix(Evirt(isym), eigen(:,nocc+1:), context1, block)
      end if
      call DeleteSharedArray(eigen,'eigen')
      deallocate(npart,tmp)
   end do
   call kfclfl(iu21)
   call NewSharedArray(gPoltddft%dipole,'dipole',(/gPoltddft%dimecc,3/))
   call NewSharedArray(gPoltddft%dipvel,'dipvel',(/gPoltddft%dimecc,3/))
   call NewSharedArray(gPoltddft%magnetic,'magnetic',(/gPoltddft%dimecc,3/))
   if (AmIOwner(gPoltddft%dipole)) gPoltddft%dipole = 0.0_KREAL
   if (AmIOwner(gPoltddft%dipvel)) gPoltddft%dipvel = 0.0_KREAL
   if (AmIOwner(gPoltddft%magnetic)) gPoltddft%magnetic = 0.0_KREAL
   call ppNodeBarrier

   call kfopfl (iu15, 'TAPE15')
   call kfopsc (iu15, 'Matrices')

   call DoOneMatrix('Dipmat_x', gPoltddft%dipole, 1, .true.)
   call DoOneMatrix('Dipmat_y', gPoltddft%dipole, 2, .true.)
   call DoOneMatrix('Dipmat_z', gPoltddft%dipole, 3, .true.)
   call DoOneMatrix('Nablamat_x', gPoltddft%dipvel, 1, .false.)
   call DoOneMatrix('Nablamat_y', gPoltddft%dipvel, 2, .false.)
   call DoOneMatrix('Nablamat_z', gPoltddft%dipvel, 3, .false.)
   call DoOneMatrix('Lmat_x', gPoltddft%magnetic, 1, .false.)
   call DoOneMatrix('Lmat_y', gPoltddft%magnetic, 2, .false.)
   call DoOneMatrix('Lmat_z', gPoltddft%magnetic, 3, .false.)

   call kfclfl(iu15)
   do isym = 1, nsymr
      call DeleteDistributedMatrix(Eocc(isym))
      call DeleteDistributedMatrix(Evirt(isym))
   enddo
   call timere('pp')

   call adf_end_blacs(context1)
   call timere ('dipole')

   contains

   subroutine DoOneMatrix(name, shared, idx, symmetric)
      character(*), intent(in)   :: name
      real(KREAL), pointer       :: shared(:,:)
      integer(KINT), intent(in)  :: idx
      logical, intent(in)        :: symmetric
!
!     Transform a matrix from ao-ao basis to virt-occ and copy selected elements to the target shared array
!
      type(DistributedMatrixType) :: virtOccDist, aoao
      real(KREAL), pointer        :: virtOcc(:,:)
      integer(KINT)               :: i, j, isymocc, nocc, isymvirt, nvirt, nrpq, idxecc

      if (symmetric) then
         call ReadSymmAOMatrixToDist(iu15, name, aoao, gDims%naos, 0, context1, block)
      else
         call ReadAntisymAOMatrixToDist(iu15, name, aoao, gDims%naos, 0, context1, block)
      end if
      idxecc = 0
      do isymocc = 1, nsymr
         nocc = nrocc(isymocc, 1)
         if (nocc == 0) cycle 
         do isymvirt = 1, nsymr
            nvirt = nrvirt(isymvirt, 1)
            irrep_: do irrep = 1, gSymmetryInfo%nrep
               if (irrepop(irrep)==0) cycle irrep_
               nrpq = gSymmetryInfo%icgseries(irrepop(irrep), gSymrespstan%isymadf2new(isymocc),&
                    & gSymrespstan%isymadf2new(isymvirt))
               if (nrpq == 0) cycle irrep_

               call CalcFromSimilarityTransform(virtOccDist, aoao, Evirt(isymvirt), Eocc(isymocc), .true., .false.)
               call NewSharedArray(virtOcc, 'DoOneMatrix:virtOcc', (/nvirt,nocc/))
               if (AmIOwner(virtOcc)) virtOcc = 0.0_KREAL
               if (IsSharedArray(virtOcc)) call ppNodeBarrier
               call ToFullMatrix(virtOccDist, virtOcc)
               call DeleteDistributedMatrix(virtOccDist)

               if (AmIOwner(shared)) then
                  do i = gPoltddft%noccmin(isymocc,isymvirt), nocc
                     do j = 1, gPoltddft%nvirtmax(i,isymocc,isymvirt)
                        idxecc = idxecc + 1
                        shared(idxecc,idx) = virtOcc(j,i)
                     enddo
                  enddo
               end if
               call DeleteSharedArray(virtOcc, 'DoOneMatrix:virtOcc')
               exit irrep_

            end do irrep_
         end do
      end do
      call DeleteDistributedMatrix(aoao)

   end subroutine
End subroutine
