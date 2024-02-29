Subroutine CreateKGrid(noccmx, nunmax, nsymr, nrocc, nrvirt, epsocc, epsun, irrepop, fragments)

   use PoltddftModule
   use ppCommonsModule   
   use Vartypes
   use ExcitSOVars
   use ExcitOpenMod
   use SymmetryInfo
   use FragmentsArrayTypeModule
   use OrbitalsCoulombInteractionModuleReduced
   use SymptrcModule
   use OrbitalsTypeModule
   use SymrespstanModule, only: gSymrespstan
   use ADFGlobalInputModule !prova max number of integrals
   use NucleiModule
   use FitsymModule
   use ADFFilesModule
   use adf_blacs
   use DistributedMatrixModule
   use SharedArraysModule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PIERPAOLO
     use fitfitint2RangeSeparatedModule
   use XCRangeSeparatedCalcDescriptor
   use XCFunctionalDescriptorModule
   use RangeSeparatedCalculation
   use XCRangeSeparatedCalcDescriptor
   use SharedArraysModule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none 

   integer(KINT), intent(in) :: noccmx, nunmax
   integer(KINT), intent(in) :: nsymr, irrepop(gSymmetryInfo%nrep)
   integer(KINT), intent(in) :: nrocc(gExcitSOVars%nsymt,gExcitOpenMod%nds)
   integer(KINT), intent(in) :: nrvirt(gExcitSOVars%nsymt,gExcitOpenMod%nds)

   real(KREAL),intent(in)    :: epsocc(noccmx,gExcitOpenMod%nds), epsun(nunmax,gExcitOpenMod%nds)
   type(FragmentsArrayType), intent(inout) :: fragments

   integer(KINT) :: inocc, invirt, i, z, z_PIER, irrep, y, delta_noccmin, fake_noccmin(nsymr,nsymr)
   integer(KINT) :: iocc, ivirt, nocc, nvirt, nrpq, min_noccmin, c, d
   integer(KINT) :: isymocc, isymvirt
   integer(KINT) :: noccmxk
   integer(KINT) :: iOrbital, jOrbital
   integer(KINT) :: index_orb(gDims%nsot,gDims%nsym)
   real(KREAL)   :: correctionenergy, hdacutoffrestart
   real(KREAL)   :: fhf, fxlda, fxgga, fclda, fcgga
   real(KREAL)   :: eps, eps_postcorrection

   integer(KINT) :: iexist64
   integer(KINT) :: iu64, iu, iu65, iu21, iu65tocc, iu65tvirt
   real(KREAL), allocatable :: epsstored(:)
   character(LCHARS) :: line !prova number of integrals
   integer(KINT) :: k = 10, kk !prova counter
   integer(KINT) :: nhdaintegral  !prova max number of integrals
   integer(KINT) :: abs_epsstored_index !prova index
   integer(KINT) :: context
   integer(KINT) :: block
   integer(KINT) :: print_it
   real(KREAL), pointer :: Q(:,:)
   real(KREAL) :: docc(1,gDims%nsfos), dvirt(1,gDims%nsfos)
   real(KREAL) :: doccQ(gDims%nsfos,1),doccQdvirt(1,1), delta_ia
   character(LCHARS) :: numberocc, numbervirt, symrep(nsymr)
   logical :: HDA_fitted, IsRangeSeparated !provaP
   type(DistributedMatrixType)   :: Q_dist
   logical :: nhdaintegral_upward = .false., nhdaintegral_overshoot = .false. !prova restart


   call timers ('CreateKgrid')

   call msg ('CREATEKGRID')

   allocate(gPoltddft%irrepdip(gSymmetryInfo%nrep))
   gPoltddft%irrepdip = 0_KINT
   gPoltddft%irrepdip(1:gSymmetryInfo%nrep) = irrepop

   allocate(gPoltddft%noccmin(nsymr,nsymr)) 
   gPoltddft%noccmin = 0
   gPoltddft%dimecc = 0
   noccmxk = 0 

   HDA_fitted = .false.
   if (gModel%lhybrid) then
      IsRangeSeparated = IsHFRangeSeparatedFunctional(gModel%xcDescriptor) !!!!!!!!!!!!!!PIERPAOLO
      call xchybpars (gModel%hybrid, fhf, fxlda, fxgga, fclda, fcgga)
      if (IsRangeSeparated) then
          fxlda = 0
          fxgga = 0
          fclda = 0
          fcgga = 0
      endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call gInput%Get ('PolTDDFT%N_HDA_integral', nhdaintegral) !prova max number of integrals
      call gInput%Get ('PolTDDFT%Print_Int_time', print_it)
      call gInput%Get ('PolTDDFT%HDA_fitted', HDA_fitted)
   end if


   iexist64 = 0
   if (gModel%lhybrid) then
      if (kfexfl('TAPE64')) iexist64 = 1
      call ppNodeBarrier
      call ppcbi(iexist64,'iexst64')
   end if

   
   if ((iexist64 == 1).and.(HDA_fitted)) then
      HDA_fitted = .false.
      call msg('HDA_fitted is not available for restart, the simulation keeps &
                going using standard numerical HDA')
      call ppNodeBarrier
   end if

   iocc = 0
   do isymocc =1, nsymr
      nocc = nrocc(isymocc, 1)
 
      if (nocc > noccmxk) noccmxk = nocc

      ivirt = 0
      do isymvirt = 1, nsymr
         nvirt = nrvirt(isymvirt, 1)

         irrep_: do irrep = 1, gSymmetryInfo%nrep
            if (irrepop(irrep)==0) cycle irrep_
            nrpq = gSymmetryInfo%icgseries(irrepop(irrep), gSymrespstan%isymadf2new(isymocc),&
                 & gSymrespstan%isymadf2new(isymvirt))

            if (nrpq == 0) cycle irrep_

            gPoltddft%noccmin(isymocc,isymvirt) = 1 

            do inocc = 1, nocc
               do invirt = 1, nvirt
                  eps = epsun(ivirt+invirt,1) - epsocc(iocc+inocc,1)
                  if (eps > gPoltddft%eVgrid) then
                     if (invirt == 1) gPoltddft%noccmin(isymocc,isymvirt) = gPoltddft%noccmin(isymocc,isymvirt) + 1 
                     exit
                  endif
                  gPoltddft%dimecc = gPoltddft%dimecc + 1               
               enddo
            enddo
            exit irrep_

         enddo irrep_

         ivirt = ivirt + nvirt
      enddo
      iocc = iocc + nocc
   enddo

   allocate (gPoltddft%kidx(gPoltddft%dimecc))
   allocate (gPoltddft%energy(gPoltddft%dimecc))
   allocate (gPoltddft%nvirtmax(noccmxk,nsymr,nsymr))
   allocate (gPoltddft%eiocc(gPoltddft%dimecc))
   allocate (gPoltddft%eivirt(gPoltddft%dimecc))

   allocate (gPoltddft%symocc(gPoltddft%dimecc))
   allocate (gPoltddft%noocc(gPoltddft%dimecc))
   allocate (gPoltddft%symvir(gPoltddft%dimecc))
   allocate (gPoltddft%novir(gPoltddft%dimecc))

   gPoltddft%nvirtmax = 0
   gPoltddft%kidx = 0
   gPoltddft%energy = 0
   gPoltddft%eiocc = 0.0_KREAL
   gPoltddft%eivirt = 0.0_KREAL

   gPoltddft%symocc = 0
   gPoltddft%noocc = 0
   gPoltddft%symvir = 0
   gPoltddft%novir = 0


   iexist64 = 0
   if (gModel%lhybrid) then
      if (kfexfl('TAPE64')) iexist64 = 1
      call ppNodeBarrier
      call ppcbi(iexist64,'iexst64')

      allocate(epsstored(noccmx*nunmax))
      epsstored = 0.0_KREAL

      if (iexist64 /= 0) then
         call ppNodeBarrier
         call ppcbi(iexist64,'iexst64')
         call kfopflpp(iu64, 'TAPE64')
         call kfclfl(iu64)

         call kfopfl(iu64, 'TAPE64')
         call kfopsc(iu64, 'diagonalexchangeB3LYP') !prova changed name section from CreateKGridB3LYP to diagonalexchangeB3LYP
         call kfread(iu64, 'epsstored', epsstored)
         call kfread(iu64, 'hdacutoff', hdacutoffrestart)
         call kfclfl(iu64)

      end if

      i = 1
      do isymocc = 1, gDims%nsym
         do iocc = 1, gSymptrc%norboc(1,isymocc)
            index_orb(iocc,isymocc) = i
            i = i + 1
         enddo
      enddo
      do isymvirt = 1, gDims%nsym
         do ivirt = gSymptrc%norboc(1,isymvirt)+1, gSymptrc%norb(isymvirt)
            index_orb(ivirt,isymvirt) = i
            i = i + 1
         enddo
      enddo

      if (.not. HDA_fitted) then !provaP

         if (iexist64 == 0 .or. (iexist64 /= 0 .and. hdacutoffrestart < gPoltddft%eVgrid)) allocate(giBlock(1))

      else

         iocc = 0
         do isymocc =1, nsymr
            nocc = nrocc(isymocc, 1)

            if (nocc == 0) cycle

            ivirt = 0
            do isymvirt = 1, nsymr
               nvirt = nrvirt(isymvirt, 1)

               irrep2_: do irrep = 1, gSymmetryInfo%nrep
                  if (irrepop(irrep)==0) cycle irrep2_
                  nrpq = gSymmetryInfo%icgseries(irrepop(irrep),gSymrespstan%isymadf2new(isymocc),&
                    & gSymrespstan%isymadf2new(isymvirt))

                  if (nrpq == 0) cycle irrep2_

                  do inocc = gPoltddft%noccmin(isymocc,isymvirt), nocc

                     gPoltddft%nvirtmax(inocc,isymocc, isymvirt) = nvirt

                     do invirt = 1, nvirt
                        eps = epsun(ivirt+invirt,1) - epsocc(iocc+inocc,1)
                        if (eps > gPoltddft%eVgrid) then

                           gPoltddft%nvirtmax(inocc,isymocc, isymvirt) = invirt - 1

                        exit
                        endif
                     enddo
                  enddo
               enddo irrep2_
               ivirt = ivirt + nvirt
            enddo
            iocc = iocc + nocc
         enddo

         call ppNodeBarrier
         call Calcdpadqi(nsymr, nrocc, irrepop) !provaP

         call ppNodeBarrier

      end if !provaP

   end if

   if ((gModel%lhybrid).and.(print_it.eq.1)) then !prova number of integrals
     write(line,*) gPoltddft%dimecc 
     line = "Total number of HDA integrals = " // trim(line)
     call msg(trim(line)) 
   end if !prova number of integrals

!INIZIO provaP
   fake_noccmin = 0 !provaP
   fake_noccmin = gPoltddft%noccmin !provaP
   do c = 1, nsymr !provaP
      do d = 1, nsymr !provaP
         if (fake_noccmin(c,d) == 0) fake_noccmin(c,d) = 1000000 !provaP
      enddo !provaP
   enddo !provaP
!FINE provaP

   kk = 0 !prova number of integrals
   if (.not. HDA_fitted) then !provaP
      gOrbitalsCoulombInteractionReduced%icounter = 0 !provaP
   else
      call GetBLACSContextAndBlocksize(gDims%nsfos, context, block)
      call kfopfl(iu65, 'TAPE65') !provaP
      call kfopsc(iu65,'Q') !provaP
      call NewDistributedMatrix(Q_dist,REAL_MATRIX_TYPE,gDims%nsfos,gDims%nsfos,context,block)
      call SetConstant(Q_dist, 0.0_KREAL)
      call LoadDistributedMatrix (Q_dist, iu65, 1, 'Q')
      call NewSharedArray(Q, 'Q', [gDims%nsfos,gDims%nsfos])
      call kfclfl(iu65) !provaP
      call ToFullMatrix (Q_dist, Q)
      call DeleteDistributedMatrix(Q_dist)

      call kfopfl (iu21, gADFFiles%main)
      call kfopvr (iu21, 'Symmetry%nsym')
      call kfread (iu21, 'symlab', symrep)
      call kfclfl (iu21)

      call kfopfl(iu65tocc, 'TAPE65t') !provaP
      call kfopfl(iu65tvirt, 'TAPE65t') !provaP

   end if
   z = 0
!write(*,*) "PIERPAOLO z index has been initialized to 0, z = ", z
   iocc = 0
   do isymocc =1, nsymr
      y = 0
      nocc = nrocc(isymocc, 1)

      if (nocc == 0) cycle
      min_noccmin = minval(fake_noccmin(isymocc,:))
      ivirt = 0
      do isymvirt = 1, nsymr
         nvirt = nrvirt(isymvirt, 1)
       
         irrep3_: do irrep = 1, gSymmetryInfo%nrep
            if (irrepop(irrep)==0) cycle irrep3_
            nrpq = gSymmetryInfo%icgseries(irrepop(irrep),gSymrespstan%isymadf2new(isymocc),&
                 & gSymrespstan%isymadf2new(isymvirt))

            if (nrpq == 0) cycle irrep3_

            do inocc = gPoltddft%noccmin(isymocc,isymvirt), nocc
               delta_noccmin = gPoltddft%noccmin(isymocc,isymvirt) - min_noccmin
               gPoltddft%nvirtmax(inocc,isymocc, isymvirt) = nvirt
               if (HDA_fitted) then !provaP
                  y = inocc - gPoltddft%noccmin(isymocc,isymvirt) + 1 + delta_noccmin  !provaP ho messo y al posto di inocc
                  numberocc = ' '
                  call csputi(numberocc, y) !(numberocc, inocc) provaP
                  call kfrdnr(iu65tocc, trim(symrep(isymocc))//'_occ%dens_fit_'//trim(numberocc), docc(1,:), gDims%nsfos, 1) !provaP
               end if !provaP
               do invirt = 1, nvirt
                  eps = epsun(ivirt+invirt,1) - epsocc(iocc+inocc,1)
                  if (eps > gPoltddft%eVgrid) then 

                     gPoltddft%nvirtmax(inocc,isymocc, isymvirt) = invirt - 1

                     exit
                  endif
                  if (HDA_fitted) then !provaP
                     numbervirt = ' '
                     call csputi (numbervirt, invirt)
                     call kfrdnr(iu65tvirt,trim(symrep(isymvirt))//'_virt%dens_fit_'//trim(numbervirt), dvirt(1,:),gDims%nsfos, 1) !provaP
                  end if
                  z = z + 1       
                  if (gModel%lhybrid) then
                     iOrbital = index_orb(inocc,isymocc)
                     jOrbital = index_orb(invirt+gSymptrc%norboc(1,isymvirt),isymvirt)
                     abs_epsstored_index = (iOrbital-1)*nunmax + jOrbital-noccmx !prova indice assoluto

                     if (epsstored(abs_epsstored_index) == 0.0_KREAL) then !prova massimo numero di integrali

                        kk = kk + 1 !prova conta integrali fatti

                        if (kk > nhdaintegral) then !prova massimo numero di integrali

                          nhdaintegral_overshoot = .true. !prova massimo numero di integrali
                          kk = kk - 1 !prova massimo numero di integrali
                          goto 10 !prova massimo numero di integrali

                        end if !prova massimo numero di integrali

                        if ((iexist64 /= 0) .and. (eps <= hdacutoffrestart)) then !prova restart numero orbitali

                           if (.not. HDA_fitted) then !provaP
                              gOrbitalsCoulombInteractionReduced%icounter = 0 !prova restart numero orbitali
                              if (.not. allocated(giBlock)) allocate(giBlock(1)) !prova restart numero orbitali
                           end if !provaP
                           nhdaintegral_upward = .true. !prova in caso di restart verso l'alto del numero di integrali
                        end if !prova restart numero orbitali
!START provaP here we use dummy matrices in order to use mmulrr, maybe it is
!possible to find a better routines.

                        if (HDA_fitted) then !provaP

                           doccQ = 0.0_KREAL
                            call mmulrr(Q,gDims%nsfos,gDims%nsfos, &
                              docc(1,:),gDims%nsfos,1,'NO',doccQ)
                           doccQdvirt = 0.0_KREAL
                            call mmulrr(dvirt(1,:),1,gDims%nsfos, &
                              doccQ(:,1),gDims%nsfos,1,'NO',doccQdvirt)
                           
                           delta_ia = doccQdvirt(1,1)*27.211
                           !eps = eps - doccQdvirt(1,1)*fhf
                           !write(*,*) "delta_ia PIERPAOLO = ", delta_ia
                           eps_postcorrection = eps - doccQdvirt(1,1)*fhf !PIERPAOLO
                           if (eps_postcorrection.lt.0.d0) then
                              eps_postcorrection = 0.d0
                           else if (eps_postcorrection.gt.gPoltddft%eVgrid) then
                              eps_postcorrection = gPoltddft%eVgrid
                           end if
                              
                           !eps = eps - doccQdvirt(1,1)*fhf
                           
                        else

!STOP provaP



                           call CalcOrbitalsCoulombInteractionReduced(fragments%activeFragment, iOrbital, jOrbital, correctionenergy)
                           eps_postcorrection = eps - correctionenergy*fhf !PIERPAOLO 
                           if (eps_postcorrection.lt.0.d0) then
                              eps_postcorrection = 0.d0
                           else if (eps_postcorrection.gt.gPoltddft%eVgrid) then
                              eps_postcorrection = gPoltddft%eVgrid
                           end if
                           !eps = eps - correctionenergy*fhf
                           
                        end if

                        epsstored(abs_epsstored_index) = eps_postcorrection !epsstored(abs_epsstored_index) = eps !PIERPAOLO epsstored(abs_epsstored_index) = eps_postcorrection
                        

                        if ((k == kk).and.(print_it.eq.1)) then !prova number of integrals
                        
                           write(line,*) kk 
                           line = "Number of HDA integrals done = " // trim(line)
                           call msg (trim(line))
                           k = k * 10
                        
                        end if !prova number of integrals

                     else
                        eps = epsstored(abs_epsstored_index) !eps = epsstored(abs_epsstored_index) !PIERPAOLO eps_postcorrection = epsstored(abs_epsstored_index)
                        eps_postcorrection = eps
                     end if
                  end if
                  
                  if (gModel%lhybrid) then ! PIER i added this if

                     if (eps <= gPoltddft%eVgrid) then !PIERPAOLO
                        do i = 1, gPoltddft%nkgrid
                           if (gPoltddft%intervals(i) <= eps_postcorrection .and. eps_postcorrection < gPoltddft%intervals(i+1)) then  !PIERPAOLO eps diventa eps_postcorrection                   
                         !if (gPoltddft%intervals(i) <= eps .and. eps < gPoltddft%intervals(i+1)) then
                              gPoltddft%kidx(z)   = i 
                              gPoltddft%energy(z) = eps_postcorrection !PIERPAOLO eps                 
                              gPoltddft%eiocc(z)  = epsocc(iocc+inocc,1)
                              gPoltddft%eivirt(z) = epsun(ivirt+invirt,1)
                              gPoltddft%symocc(z) = isymocc
                              gPoltddft%noocc(z)  = inocc
                              gPoltddft%symvir(z) = isymvirt
                              gPoltddft%novir(z)  = invirt + nrocc(isymvirt,1)
                           endif
                        enddo
                     else if (eps > gPoltddft%eVgrid) then
                        z = z-1
                  !  write(iuout,*) "PIERPAOLO eps > gPoltddft%eVgrid con z-1 =", z
                     endif   !PIERPAOLO
                  else
                     do i = 1, gPoltddft%nkgrid
                        if (gPoltddft%intervals(i) <= eps .and. eps < gPoltddft%intervals(i+1)) then

                           gPoltddft%kidx(z)   = i
                           gPoltddft%energy(z) = eps
                           gPoltddft%eiocc(z)  = epsocc(iocc+inocc,1)
                           gPoltddft%eivirt(z) = epsun(ivirt+invirt,1)

                           gPoltddft%symocc(z) = isymocc
                           gPoltddft%noocc(z)  = inocc
                           gPoltddft%symvir(z) = isymvirt
                           gPoltddft%novir(z)  = invirt + nrocc(isymvirt,1)
 
                        endif
                     enddo
                  endif
               enddo
            enddo
            exit irrep3_
         enddo irrep3_
         ivirt = ivirt + nvirt
      enddo
      iocc = iocc + nocc
   enddo

!write(iuout,*) "PIERPAOLO gPoltddft%kidx", gPoltddft%kidx
!write(iuout,*) "PIERPAOLO gPoltddft%energy", gPoltddft%energy
!write(iuout,*) "PIERPAOLO gPoltddft%eiocc", gPoltddft%eiocc
!write(iuout,*) "PIERPAOLO gPoltddft%eivirt", gPoltddft%eivirt
!write(iuout,*) "PIERPAOLO gPoltddft%symocc", gPoltddft%symocc
!write(iuout,*) "PIERPAOLO gPoltddft%noocc", gPoltddft%noocc
!write(iuout,*) "PIERPAOLO gPoltddft%symvir", gPoltddft%symvir
!write(iuout,*) "PIERPAOLO gPoltddft%novir", gPoltddft%novir
!write(iuout,*) "PIERPAOLO gPoltddft%nvirtmax", gPoltddft%nvirtmax
!write(iuout,*) "PIERPAOLO dimecc", gPoltddft%dimecc
!write(iuout,*) "PIERPAOLO gPoltddft%noccmin", gPoltddft%noccmin
!write(iuout,*) "PIERPAOLO gPoltddft%irrepdip", gPoltddft%irrepdip
!write(iuout,*) "PIERPAOLO gPoltddft%nkgrid", gPoltddft%nkgrid
!write(iuout,*) "PIERPAOLO gPoltddft%intervals", gPoltddft%intervals



!write(*,*) "PIERPAOLO array element, ", z_PIER," of gPoltddft%symocc is:", gPoltddft%symocc(z_PIER)
!write(*,*) "PIERPAOLO array element, ", z_PIER," of gPoltddft%symvir is:", gPoltddft%symvir(z_PIER)
!write(*,*) "PIERPAOLO number of eccitation is gPoltddft%dimecc =", gPoltddft%dimecc

!call stopit("PIERPAOLO")

   if (HDA_fitted) then
   call kfclfl(iu65tocc)
   call kfclfl(iu65tvirt)
   endif

   allocate (gPoltddft%kdim(gPoltddft%nkgrid))
   gPoltddft%kdim = 0
   do i = 1, gPoltddft%dimecc
      if (gPoltddft%kidx(i) > 0) &
         gPoltddft%kdim(gPoltddft%kidx(i)) = gPoltddft%kdim(gPoltddft%kidx(i))+1
   enddo
   10 continue !prova max number of integrals
   if (gModel%lhybrid) then
      if (iexist64 /= 0 .and. ((gPoltddft%eVgrid > hdacutoffrestart) .or. nhdaintegral_upward)) then !prova in caso di restart verso l'alto del numero di integrali)
         if (.not. HDA_fitted) deallocate(giBlock) !provaP
         call kfopfl(iu64, 'TAPE64')
         call kfdlsc(iu64, 'diagonalexchangeB3LYP') !prova changed name section from CreateKGridB3LYP to diagonalexchangeB3LYP
         call kfcrsc(iu64, 'diagonalexchangeB3LYP') !prova changed name section from CreateKGridB3LYP to diagonalexchangeB3LYP
         call kfwrite(iu64, 'hdacutoff', gPoltddft%eVgrid)
         call kfwrnr(iu64, 'epsstored', epsstored, noccmx*nunmax, 1)
         call kfclfl(iu64)
      else if (iexist64 == 0) then
         if (.not. HDA_fitted) deallocate(giBlock) !provaP
         call kfcrfl(iu64, 'TAPE64')
         call kfcrsc(iu64, 'diagonalexchangeB3LYP') !prova changed name section from CreateKGridB3LYP to diagonalexchangeB3LYP
         call kfwrite(iu64, 'hdacutoff', gPoltddft%eVgrid)
         call kfwrnr(iu64, 'epsstored', epsstored, noccmx*nunmax, 1)
         call kfclfl(iu64)
      end if
      if (HDA_fitted) call DeleteSharedArray (Q, 'Q') !provaP pulizia matrice Q
      deallocate(epsstored)
      if (gFitsym%lnosft) then
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
   end if

   call timere ('CreateKgrid')

   if ((gModel%lhybrid).and.(print_it.eq.1)) then !prova max number of integrals
     write(line,*) kk 
     line = "Done " // trim(line) // " HDA integrals and written TAPE64"
     call msg (trim(line)) !prova max number of integrals
   end if

   if (nhdaintegral_overshoot) call stopit ('CreateKGrid: not all HDA integrals are done, you can restart from TAPE64') !prova max number of integrals

write(iuout,*) "end of CreateKGrid PIER"

End subroutine CreateKGrid 
