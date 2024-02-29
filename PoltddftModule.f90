module PoltddftModule

   use Vartypes 
   use ADFFilesModule
   implicit none
   
   type PoltddftType
      real(KREAL)   :: eVgrid, wi, lambda
      real(KREAL)   :: cutoff
      integer(KINT) :: nkgrid, npoints
      integer(KINT) :: dimecc ! Number of occ-virt MO pairs
      logical       :: dv
      integer(KINT), allocatable :: irrepdip(:)

      logical, allocatable       :: skipnb(:)
      integer(KINT), allocatable :: noccmin(:,:)
      integer(KINT), allocatable :: nvirtmax(:,:,:)
      integer(KINT), allocatable :: kidx(:)
      integer(KINT), allocatable :: kdim(:)
      integer(KINT), allocatable :: symocc(:)    ! (dimecc)
      integer(KINT), allocatable :: noocc(:)     ! (dimecc)
      integer(KINT), allocatable :: symvir(:)    ! (dimecc)
      integer(KINT), allocatable :: novir(:)     ! (dimecc)
      real(KREAL), allocatable   :: energy(:)
      real(KREAL), allocatable   :: intervals(:) ! (nkgrid)   energy at grid points???
      real(KREAL), allocatable   :: alpr(:,:)    ! (npoints,3)  
      real(KREAL), allocatable   :: alpi(:,:)    ! (npoints,3)  
      real(KREAL), allocatable   :: betr(:,:)    ! (npoints,3)  
      real(KREAL), allocatable   :: beti(:,:)    ! (npoints,3)  
      real(KREAL), allocatable   :: wr(:)        ! (npoints)  
      real(KREAL), allocatable   :: eiocc(:)     ! (dimecc)
      real(KREAL), allocatable   :: eivirt(:)    ! (dimecc)
      real(KREAL), allocatable   :: ianalysis(:) ! (dimecc)
      real(KREAL), allocatable   :: ranalysis(:) ! (dimecc)
      real(KREAL), pointer       :: dipole(:,:)   ! (dimecc, 3)  
      real(KREAL), pointer       :: dipvel(:,:)   ! (dimecc, 3)  
      real(KREAL), pointer       :: magnetic(:,:) ! (dimecc, 3)  
      real(KREAL), pointer       :: Vx(:,:) ! (nsfos,nkgrid)
      real(KREAL), pointer       :: Vy(:,:) ! (nsfos,nkgrid)
      real(KREAL), pointer       :: Vz(:,:) ! (nsfos,nkgrid)
      real(KREAL), pointer       :: Vmx(:,:) ! (nsfos,nkgrid)
      real(KREAL), pointer       :: Vmy(:,:) ! (nsfos,nkgrid)
      real(KREAL), pointer       :: Vmz(:,:) ! (nsfos,nkgrid)
   end type PoltddftType

   type(PoltddftType), save :: gPoltddft
   
   contains
   subroutine PrintandWriteKFPolTDDFT
      use PhyconModule
      use KF
      use SymmetryInfo

      real(KREAL), parameter     :: two = 2.0_KREAL, three = 3.0_KREAL

      real(KREAL), allocatable   :: polarizability(:)
      real(KREAL), allocatable   :: absorption(:)
      real(KREAL), allocatable   :: cdspectrum(:)

      integer(KINT)  :: i, iu, j, k1, k2

      real(KREAL) :: eVha
      real(KREAL), external :: convrs
      integer(KINT) :: lvaluedipole = 1
      real(KREAL) :: tensorops(3,3)
      integer(KINT) :: irreps(3), ndimirreps(3), multiplicity(3), irrepcount
      real(KREAL), parameter :: k_deltaepsilon = 0.0138607570264085_KREAL !prova [16*(PI**2)*N / 3*(h/PI)*c*ln(10)*(10**3)] / PI
      real(KREAL), parameter :: k_epsilon = 9066.63970956273_KREAL !prova {[2*(PI**2)*N / (h/PI)*c*ln(10)*(10**3)] * [eVha*(Bohr*Statcoulomb)**2]} / PI
      real(KREAL), allocatable   :: epsilon(:), deltaepsilon(:) !prova epsilon and deltaepsilon
      

      write(iuout, 9000)

      call IrredTensorOps (lvaluedipole, tensorops, irreps, ndimirreps, multiplicity, irrepcount)

      write(iuout,*)
      write(iuout,*) "Allowed symmetries:"
      write(iuout,*)

      k1 = 0
      k2 = 0

      do i = 1, 3
        if (irreps(i) /= 0) then
          write(iuout,*) gSymmetryInfo%symlab(irreps(i))
          k1 = k1 + 1
        end if
      end do

      do i = 1, 3
        if (irreps(i) == 0) cycle
        do j = 1, gSymmetryInfo%nrep
          if (gPoltddft%irrepdip(j) == irreps(i)) k2 = k2 + 1
        end do
      end do

      deallocate(gPoltddft%irrepdip)

      eVha = convrs('ENERGY','HARTREE,EV')

      write(iuout,'(1x,A)') '',&
         '=============================================================================================================================',&
         'Pol-TDDFT results'
      write(iuout,'(1x,A)') &
         '-----------------------------------------------------------------------------------------------------------------------------',&
         '   E(Hartree)         E(eV)       A(real,X)      A(imag,X)       A(real,Y)       A(imag,Y)       A(real,Z)       A(imag,Z)',&
         '-----------------------------------------------------------------------------------------------------------------------------'
      do i = 1, gPoltddft%npoints
         write(iuout,3000) gPoltddft%wr(i), gPoltddft%wr(i)*eVha, &
               gPoltddft%alpr(i,1), gPoltddft%alpi(i,1)*gPoltddft%wr(i)*2/3*gPoltddft%wi, &
               gPoltddft%alpr(i,2), gPoltddft%alpi(i,2)*gPoltddft%wr(i)*2/3*gPoltddft%wi, &
               gPoltddft%alpr(i,3), gPoltddft%alpi(i,3)*gPoltddft%wr(i)*2/3*gPoltddft%wi
      enddo

      write(iuout,'(1x,A)') &
           '-----------------------------------------------------------------------------------------------------------------------------',&
           '   E(Hartree)         E(eV)       B(imag,X)      B(imag,Y)       B(imag,Z)',&
           '-----------------------------------------------------------------------------'
      do i = 1, gPoltddft%npoints
         write(iuout,3000) gPoltddft%wr(i), gPoltddft%wr(i)*eVha, gPoltddft%beti(i,1)*gPoltddft%wi*rotau2cgs, &
                                                                  gPoltddft%beti(i,2)*gPoltddft%wi*rotau2cgs, &
                                                                  gPoltddft%beti(i,3)*gPoltddft%wi*rotau2cgs
      enddo

      if (k1 == k2) then

         allocate (polarizability(gPoltddft%npoints))
         allocate (absorption(gPoltddft%npoints))
         allocate (cdspectrum(gPoltddft%npoints))
         allocate(epsilon(gPoltddft%npoints), deltaepsilon(gPoltddft%npoints)) !prova epsilon and deltaepsilon

         do i = 1, gPoltddft%npoints
            polarizability(i) = sum(gPoltddft%alpr(i,:))/three
            absorption(i) = sum(gPoltddft%alpi(i,:))*gPoltddft%wr(i)*gPoltddft%wi*two/three
            cdspectrum(i) = sum(gPoltddft%beti(i,:))*gPoltddft%wi*rotau2cgs
            deltaepsilon(i) = k_deltaepsilon*gPoltddft%wr(i)*cdspectrum(i) / gPoltddft%wi !prova epsilon and deltaepsilon (approximation)
         enddo

         epsilon = (k_epsilon*absorption) / (gPoltddft%wi*eVha) !prova epsilon and deltaepsilon

         write(iuout,'(1x,A)') &
              '---------------------------------------------------------------------------------------------------------------', &
              '   E(Hartree)         E(eV)    Polariz.(real)    Absorption      CD spectrum      Epsilon      Delta_Epsilon',  &
              '---------------------------------------------------------------------------------------------------------------'
         do i = 1, gPoltddft%npoints
            write(iuout,3000) gPoltddft%wr(i), gPoltddft%wr(i)*eVha, polarizability(i), absorption(i), cdspectrum(i), epsilon(i), deltaepsilon(i) !prova epsilon and deltaepsilon
         enddo

         call kfopfl (iu, gADFFiles%main)
         call kfopsc (iu, 'POLTDDFT')
         call kfwrite (iu, 'Polarizability (real)', polarizability)
         call kfwrite (iu, 'Absorption', absorption)
         call kfwrite (iu, 'CD Spectrum', cdspectrum)
         call kfwrite (iu, 'Epsilon', epsilon) !prova epsilon and deltaepsilon
         call kfwrite (iu, 'Deltaepsilon', deltaepsilon) !prova epsilon and deltaepsilon
         call kfclfl (iu)
         deallocate (polarizability)
         deallocate (absorption)
         deallocate (cdspectrum)
         deallocate(epsilon,deltaepsilon) !prova epsilon and deltaepsilon

      else
         write(iuout,'(1x,A)') &
              '-----------------------------------------------------------------------------', &
              '   E(Hartree)         E(eV)    Polariz.(real)    Absorption      CD spectrum',  &
              '-----------------------------------------------------------------------------'
         write(iuout,*) 
         write(iuout,*) "PolTDDFT: Not printed because not all active dipole representation have been requested."
      end if

 3000 format (2F15.8,6E16.8)

 9000 format(//'****************************************************************************'/         &
             '*                                                                          *'/         &
             '*             Start POLTDDFT program part                                  *'/         &
             '*                                                                          *'/         &
             '* Relevant references are:                                                 *'/         &
             '*   O. Baseggio, G. Fronzoni, M. Stener,                                   *'/         &
             '*      The Journal of Chemical Physics, 143, 024106 (2015)                 *'/         &
             '*                                                                          *'/         &
             '*   O. Baseggio, M. De Vetta, G. Fronzoni, M. Stener, A. Fortunelli,       *'/         & 
             '*      International Journal of Quantum Chemistry, 116 (2016) 1603 - 1611  *'/         &
             '*                                                                          *'/         &
             '*   O. Baseggio, D. Toffoli, G. Fronzoni, M. Stener,                       *'/         &
             '*      L. Sementa, A. Fortunelli,                                          *'/         & 
             '*      The Journal of Physical Chemistry C, 2016, 120 (42), pp 24335-24345 *'/         &
             '*                                                                          *'/         &
             '*                                                                          *'/         &
             '****************************************************************************'//)


   end subroutine      

   subroutine CreateKFPolTDDFT
      use KF

      integer(KINT)  :: iu

      call kfopfl (iu, gADFFiles%main)
      if (kfexsc(iu, 'POLTDDFT')) call kfdlsc(iu, 'POLTDDFT')
      call kfcrsc (iu, 'POLTDDFT')
      call kfwrite (iu, 'number_of_frequencies', gPoltddft%npoints)
      call kfwrite (iu, 'frequencies', gPoltddft%wr(1:gPoltddft%npoints))
      call kfwrite (iu, 'im_freq', gPoltddft%wi)
      call kfwrite (iu, 'number_of_frequencies', gPoltddft%npoints)
      call kfwrite (iu, 'noccvirtpairs', gPoltddft%dimecc)

      call kfwrite (iu, 'symocc', gPoltddft%symocc)
      call kfwrite (iu, 'nocc', gPoltddft%noocc)
      call kfwrite (iu, 'symvir', gPoltddft%symvir)
      call kfwrite (iu, 'nvir', gPoltddft%novir)

      call kfwrite (iu, 'eiocc', gPoltddft%eiocc)
      call kfwrite (iu, 'eivirt', gPoltddft%eivirt)
      call kfclfl (iu)

   end subroutine 

end module PoltddftModule
