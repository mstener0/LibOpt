Subroutine ReorderKGrid(iexist)

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
   use PairsModule
   use MasterModule

   implicit none 

   integer(KINT), intent(in)  :: iexist

   integer(KINT)     :: i
   integer(KINT)     :: iu63, iu63b
   character(LCHARS) :: reorderstr, istr

   real(KREAL), parameter :: TOL = 1E-20

   integer(KINT), allocatable :: reorder(:)
   integer(KINT), allocatable :: tmpi(:)

   real(KREAL), allocatable   :: tmp(:)
   real(KREAL), allocatable   :: tempMat1(:),tempMat2(:)   

   call timers ('reorder')

   allocate(tmp(gPoltddft%dimecc))
   allocate(tmpi(gPoltddft%dimecc))
   allocate (reorder(gPoltddft%dimecc))

   tmp(1:gPoltddft%dimecc) = gPoltddft%energy(1:gPoltddft%dimecc)
   call sorxr (gPoltddft%dimecc, tmp, TOL, reorder)

   gPoltddft%energy(1:gPoltddft%dimecc) = tmp(reorder(1:gPoltddft%dimecc))
   tmp(1:gPoltddft%dimecc) = gPoltddft%eiocc(1:gPoltddft%dimecc)
   gPoltddft%eiocc(1:gPoltddft%dimecc) = tmp(reorder(1:gPoltddft%dimecc))
   tmp(1:gPoltddft%dimecc) = gPoltddft%eivirt(1:gPoltddft%dimecc)
   gPoltddft%eivirt(1:gPoltddft%dimecc) = tmp(reorder(1:gPoltddft%dimecc))
   tmpi(1:gPoltddft%dimecc) = gPoltddft%kidx(1:gPoltddft%dimecc)
   gPoltddft%kidx(1:gPoltddft%dimecc) = tmpi(reorder(1:gPoltddft%dimecc))

   tmpi(1:gPoltddft%dimecc) = gPoltddft%symocc(1:gPoltddft%dimecc)
   gPoltddft%symocc(1:gPoltddft%dimecc) = tmpi(reorder(1:gPoltddft%dimecc))
   tmpi(1:gPoltddft%dimecc) = gPoltddft%noocc(1:gPoltddft%dimecc)
   gPoltddft%noocc(1:gPoltddft%dimecc) = tmpi(reorder(1:gPoltddft%dimecc))
   tmpi(1:gPoltddft%dimecc) = gPoltddft%symvir(1:gPoltddft%dimecc)
   gPoltddft%symvir(1:gPoltddft%dimecc) = tmpi(reorder(1:gPoltddft%dimecc))
   tmpi(1:gPoltddft%dimecc) = gPoltddft%novir(1:gPoltddft%dimecc)
   gPoltddft%novir(1:gPoltddft%dimecc) = tmpi(reorder(1:gPoltddft%dimecc))

   deallocate(tmpi)

   if (AmIOwner(gPoltddft%dipole)) then
      tmp(1:gPoltddft%dimecc) = gPoltddft%dipole(1:gPoltddft%dimecc,1)
      gPoltddft%dipole(1:gPoltddft%dimecc,1) = tmp(reorder(1:gPoltddft%dimecc))
      tmp(1:gPoltddft%dimecc) = gPoltddft%dipole(1:gPoltddft%dimecc,2)
      gPoltddft%dipole(1:gPoltddft%dimecc,2) = tmp(reorder(1:gPoltddft%dimecc))
      tmp(1:gPoltddft%dimecc) = gPoltddft%dipole(1:gPoltddft%dimecc,3)
      gPoltddft%dipole(1:gPoltddft%dimecc,3) = tmp(reorder(1:gPoltddft%dimecc))
   end if
   if (AmIOwner(gPoltddft%dipvel)) then
      tmp(1:gPoltddft%dimecc) = gPoltddft%dipvel(1:gPoltddft%dimecc,1)
      gPoltddft%dipvel(1:gPoltddft%dimecc,1) = tmp(reorder(1:gPoltddft%dimecc))
      tmp(1:gPoltddft%dimecc) = gPoltddft%dipvel(1:gPoltddft%dimecc,2)
      gPoltddft%dipvel(1:gPoltddft%dimecc,2) = tmp(reorder(1:gPoltddft%dimecc))
      tmp(1:gPoltddft%dimecc) = gPoltddft%dipvel(1:gPoltddft%dimecc,3)
      gPoltddft%dipvel(1:gPoltddft%dimecc,3) = tmp(reorder(1:gPoltddft%dimecc))
   end if
   if (AmIOwner(gPoltddft%magnetic)) then
      tmp(1:gPoltddft%dimecc) = gPoltddft%magnetic(1:gPoltddft%dimecc,1)
      gPoltddft%magnetic(1:gPoltddft%dimecc,1) = tmp(reorder(1:gPoltddft%dimecc))
      tmp(1:gPoltddft%dimecc) = gPoltddft%magnetic(1:gPoltddft%dimecc,2)
      gPoltddft%magnetic(1:gPoltddft%dimecc,2) = tmp(reorder(1:gPoltddft%dimecc))
      tmp(1:gPoltddft%dimecc) = gPoltddft%magnetic(1:gPoltddft%dimecc,3)
      gPoltddft%magnetic(1:gPoltddft%dimecc,3) = tmp(reorder(1:gPoltddft%dimecc))
   end if
   call ppNodeBarrier

   deallocate(tmp)

   if (iexist == 0) then

      call kfopfl(iu63, 'TAPE63t')
      call kfopsc(iu63, 'MatrixAkt')
      call kfcrfl(iu63b, 'TAPE63')
      call kfcrsc(iu63b, 'MatrixAk')

      allocate(tempMat1(gDims%nsfos))
      allocate(tempMat2(gDims%nsfos))

      do i = 1, gPoltddft%dimecc
         reorderstr = ' '
         call csputi (reorderstr, reorder(i))
         call kfread (iu63, 'Ak'//trim(reorderstr), tempMat1)
         istr = ' '
         call csputi (istr, i)
         call kfwrnr(iu63b, 'Ak'//trim(istr), tempMat1, gDims%nsfos, 1)
      enddo
      call kfclfl(iu63)
      call kfdlfl('TAPE63t')
      call kfclfl(iu63b)

      deallocate(tempMat1)
      deallocate(tempMat2)

   endif

   deallocate (reorder)

   call timere ('reorder')

End subroutine ReorderKGrid 
