Subroutine CalcVxyz()

   use Vartypes
   use KF
   use PoltddftModule
   use ppCommonsModule   
   use ExcitSOVars
   use ExcitOpenMod
   use SymmetryInfo
   use DimensionsModule
   use SharedArraysModule
   use SharedArraysUtilModule
   use PairsModule
   use MasterModule

   implicit none 

   integer(KINT)     :: i, ki, kd, ilow, iup, k, n
   integer(KINT)     :: iu63
   character(LCHARS) :: kstr
  
   real(KREAL), pointer     :: dipole(:,:), magnetic(:,:)
   real(KREAL), allocatable :: Ak(:,:)

   call timers ('Vxyz')

   call NewSharedArray(gPoltddft%Vx, 'Vx', [gDims%nsfos, gPoltddft%nkgrid])
   call NewSharedArray(gPoltddft%Vy, 'Vy', [gDims%nsfos, gPoltddft%nkgrid])
   call NewSharedArray(gPoltddft%Vz, 'Vz', [gDims%nsfos, gPoltddft%nkgrid])
   call NewSharedArray(gPoltddft%Vmx, 'Vmx', [gDims%nsfos,gPoltddft%nkgrid])
   call NewSharedArray(gPoltddft%Vmy, 'Vmy', [gDims%nsfos, gPoltddft%nkgrid])
   call NewSharedArray(gPoltddft%Vmz, 'Vmz', [gDims%nsfos, gPoltddft%nkgrid])

   call SetSharedArrayToZero(gPoltddft%Vx)
   call SetSharedArrayToZero(gPoltddft%Vy)
   call SetSharedArrayToZero(gPoltddft%Vz)
   call SetSharedArrayToZero(gPoltddft%Vmx)
   call SetSharedArrayToZero(gPoltddft%Vmy)
   call SetSharedArrayToZero(gPoltddft%Vmz)

   if (gPoltddft%dv) then
      dipole => gPoltddft%dipvel
   else
      dipole => gPoltddft%dipole
   endif
   magnetic => gPoltddft%magnetic

   call kfopfl(iu63, 'TAPE63')
   call kfopsc(iu63, 'MatrixAk')

   call timers ('pp')
   ilow = 1
   iup = gPoltddft%nkgrid
   if (IsSharedArray(gPoltddft%Vx)) &
      call ppNodeDistrib(gPoltddft%nkgrid, ilow, iup)

   n  = 0
   ki = 1
   do i = 1, iup !gPoltddft%nkgrid
      kd = gPoltddft%kdim(i)
      if (kd == 0) cycle

      if (ilow <= i .and. i <= iup) then
         allocate(Ak(gDims%nsfos,kd))
         do k = n+1, n+kd
           kstr = ' '
           call csputi (kstr, k)
           call kfread (iu63, 'Ak'//trim(kstr), Ak(1:gDims%nsfos, k-n))
         end do

         gPoltddft%Vx(1:gDims%nsfos,i) = matmul(Ak(1:gDims%nsfos,1:kd),dipole(ki:ki+kd-1,1))
         gPoltddft%Vy(1:gDims%nsfos,i) = matmul(Ak(1:gDims%nsfos,1:kd),dipole(ki:ki+kd-1,2))
         gPoltddft%Vz(1:gDims%nsfos,i) = matmul(Ak(1:gDims%nsfos,1:kd),dipole(ki:ki+kd-1,3))
         gPoltddft%Vmx(1:gDims%nsfos,i) = matmul(Ak(1:gDims%nsfos,1:kd),magnetic(ki:ki+kd-1,1))
         gPoltddft%Vmy(1:gDims%nsfos,i) = matmul(Ak(1:gDims%nsfos,1:kd),magnetic(ki:ki+kd-1,2))
         gPoltddft%Vmz(1:gDims%nsfos,i) = matmul(Ak(1:gDims%nsfos,1:kd),magnetic(ki:ki+kd-1,3))
         deallocate(Ak)
      end if

      n  = n + kd
      ki = ki + gPoltddft%kdim(i)
   enddo
   call ppNodeBarrier
   call timere ('pp')

   call kfclfl(iu63)

   call timere ('Vxyz')

End subroutine CalcVxyz
