Subroutine CalcAk(nsymr, nrocc, nrvirt, irrepop, LEND63t)

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
   use ADFFilesModule
   use SymrespstanModule, only: gSymrespstan
   use PairsModule
   use MasterModule
   use ADFGlobalInputModule !prova max number of submatrices
   use RestartdataModule

   use AssertionsModule   

   implicit none 

   logical, intent(inout)     :: LEND63t
   integer(KINT), intent(in)  :: nsymr, irrepop(gSymmetryInfo%nrep)
   integer(KINT), intent(in)  :: nrocc(gExcitSOVars%nsymt,gExcitOpenMod%nds),&
                                 nrvirt(gExcitSOVars%nsymt,gExcitOpenMod%nds)

   integer(KINT) :: i, j, l, irrep, k, n
   integer(KINT) :: jvirt, jocc
   integer(KINT) :: jatyp, ja
   integer(KINT) :: nocc, nvirt, virt1, nrpq 
   integer(KINT) :: isymocc, isymvirt
   integer(KINT) :: virtatom, occatom, maxatom, minatom
   integer(KINT) :: virtbas, occbas, maxbas, minbas
   integer(KINT) :: lalo, lahi, maa
   integer(KINT) :: basa, basb, cplbas
   integer(KINT) :: ipair
   integer(KINT) :: iu21, iu62, iu63t !prova change from t63 to t63t
   integer(KINT) :: atombas(gDims%nnuc+1)
   integer(KINT) :: nbas(nsymr), nbasmax, nsym, nmo 

   character(LCHARS) :: symrep(nsymr)

   integer(KINT), allocatable :: npart(:,:), invbas(:,:)
   integer(KINT), allocatable :: ibasa(:), ibasb(:)
   integer(KINT), allocatable :: skip(:), idxpf(:)

   real(KREAL), allocatable   :: pairfitkf(:)
   real(KREAL), allocatable   :: pfcol(:)

   character(LCHARS)          :: pairlabel

   real(KREAL), allocatable   :: Ak(:,:)

   real(KREAL), allocatable   :: Aktmp(:)
   character(LCHARS)          :: kstr

   real(KREAL), pointer       :: eigenocc(:,:), eigenvirt(:,:)
   real(KREAL), pointer       :: pairfit(:,:)
   character(LCHARS) :: line, pline !prova number of submatrices
   integer(KINT) :: iexist63t, threshold, n_submatrices_restart, kk, iexist21 !prova restart
   integer(KINT) :: counter, print_it !prova counter

!prova delete iocc, ivirt and idxecc because useless
write(iuout,*) "start of CalckAk PIER"

   call timers ('Ak')

   call gInput%Get ('PolTDDFT%N_SubMatricesAk', threshold) !prova threshold
   call gInput%Get ('PolTDDFT%Print_Int_time', print_it)
   
   iexist63t = 0
   iexist21 = 0
   if (kfexfl('TAPE63t')) iexist63t = 1
   if (kfexfl('TAPE21')) iexist21 = 1

   if ((iexist63t.eq.1).and.(iexist21.eq.1)) then
      call kfopfl(iu21, gRestartdata%rsfile)
   end if
   if ((iexist63t.eq.1).and.(iexist21.eq.0)) then
      call kfopfl(iu21, gADFFiles%main)
      call msg('WARNING: you should restart from TAPE21 generated from the same &
                simulation which generated TAPE63t in order to guarantee orbital phase &
                consistency')
   end if
   if (iexist63t.eq.0) then
      call kfopfl(iu21, gADFFiles%main)
   end if
 
   call kfopvr (iu21, 'Symmetry%nsym')
   call kfread (iu21, 'Symmetry%nsym', nsym)      
   call kfread (iu21, 'symlab', symrep)
   call kfread (iu21, 'nfcn', nbas)

   nbasmax = 1
   do i = 1, nsymr
      if (nbas(i)>nbasmax) nbasmax = nbas(i)
   enddo

   allocate (npart(nbasmax,nsymr))
   allocate (invbas(nbasmax,nsymr))
   npart = 0 
   invbas = 0

   do i = 1, nsymr
      call kfopsc(iu21, symrep(i))
      call kfread(iu21, 'npart', npart(1:nbas(i),i))
   enddo

   i = 0 
   atombas(1) = 1
   do jatyp = 1, gDims%ntyp
      lalo = gTypesc%nbptr(jatyp)
      lahi = gTypesc%nbptr(jatyp+1) - 1
      maa  = lahi - lalo + 1
      do ja = gTypesc%nqptr(jatyp), gTypesc%nqptr(jatyp+1) - 1
         i = i + 1
         atombas(i+1) = atombas(i) + maa
      enddo
   enddo

   do jatyp = 1, nsymr
      do i = 1, nbas(jatyp)
         do j = 1, gDims%nnuc    
            if (npart(i,jatyp)>=atombas(j) .and. npart(i,jatyp)<atombas(j+1))then
               invbas(i,jatyp) = j
               exit
            endif
         enddo
      enddo
   enddo

   if (print_it.eq.1) then
      write(line,*) gPoltddft%dimecc !prova max number of submatrices
      line = "Total number of Ak submatrices = " // trim(line) !prova max number of submatrices
      call msg(trim(line)) !prova max number of submatrices
   endif

   ipair = 0 

   call timers('pp')

   iexist63t = 0 !prova restart
   if (kfexfl('TAPE63t')) iexist63t = 1 !prova restart
   call ppNodeBarrier !prova restart
   call ppcbi(iexist63t,'iexst63t') !prova restart

   n_submatrices_restart = 0

   if (iexist63t /= 0) then !prova restart

      call kfopflpp(iu63t, 'TAPE63t') 
      call kfclfl(iu63t)

      call kfopfl(iu63t, 'TAPE63t')
      call kfopsc(iu63t, 'MatrixAkt') 
      call kfread(iu63t, 'MatrixAkt%n_submatrices_restart', n_submatrices_restart)
      if (print_it.eq.1) then
         write(pline,*) n_submatrices_restart !prova max number of submatrices
         pline = "Number of submatrices Ak already stored in TAPE63t = " // trim(pline) !prova max number of submatrices
         call msg(trim(pline)) !prova max number of submatrices
      end if
   else !prova restart

     call kfcrfl(iu63t, 'TAPE63t') !prova change from t63 to t63t
     call kfcrsc(iu63t, 'MatrixAkt') !prova change from t63 to t63t

   end if 

   k = n_submatrices_restart !prova restart
   threshold = threshold + n_submatrices_restart !prova restart
   allocate (pfcol(gDims%nsfos))
   pfcol = 0.0_KREAL  
   kk = 0 !prova restart 
   counter = 0 !prova restart
   call kfopfl (iu62, 'TAPE62') !prova moved
   do isymocc = 1, nsymr
      nocc = nrocc(isymocc, 1)
      if (nocc == 0) cycle !prova logic substitution

      call kfopsc(iu21, symrep(isymocc))
      call NewSharedArray(eigenocc,'eigenocc',(/nbas(isymocc),nocc/))
      if (AmIOwner(eigenocc)) call kfread(iu21, 'Eigen-Bas_A', eigenocc)
      if (IsSharedArray(eigenocc)) call ppNodeBarrier
      do isymvirt = 1, nsymr
         nvirt = nrvirt(isymvirt, 1)
         virt1 = nrocc(isymvirt, 1)
         irrep_: do irrep = 1, gSymmetryInfo%nrep 
            if (irrepop(irrep)==0) cycle irrep_
            nrpq = gSymmetryInfo%icgseries(irrepop(irrep),gSymrespstan%isymadf2new(isymocc),&
                 & gSymrespstan%isymadf2new(isymvirt))

            if (nrpq == 0) cycle irrep_
            call kfopsc(iu21, symrep(isymvirt))
            call kfread(iu21, 'nmo_A', nmo)        
            call NewSharedArray(eigenvirt,'eigenvirt',(/nbas(isymvirt),nmo/))
            if (AmIOwner(eigenvirt)) call kfread(iu21, 'Eigen-Bas_A', eigenvirt)
            if (IsSharedArray(eigenvirt)) call ppNodeBarrier
            call NewSharedArray(pairfit,'pairfit',(/gDims%nsfos,nbas(isymvirt)/))
            do i = gPoltddft%noccmin(isymocc,isymvirt), nocc
               if (AmIOwner(pairfit)) pairfit = 0.0_KREAL
               if (IsSharedArray(pairfit)) call ppNodeBarrier
               if ((kk + gPoltddft%nvirtmax(i,isymocc,isymvirt)) <= n_submatrices_restart) then !prova restart

                 kk = kk + gPoltddft%nvirtmax(i,isymocc,isymvirt) !prova restart
                 cycle !prova restart

               end if !prova restart
               if ((k + gPoltddft%nvirtmax(i,isymocc,isymvirt)) > threshold) then !prova restart
                 call kfwrite(iu63t, 'n_submatrices_restart', k) 
                 call kfclfl(iu63t) 
                 call kfclfl(iu62)
                 call kfdlfl('TAPE62')
                 call kfclfl(iu21)
                 if (IsSharedArray(pairfit)) call ppNodeBarrier
                 call DeleteSharedArray(pairfit,'pairfit')
                 call DeleteSharedArray(eigenvirt,'eigenvirt')
                 call DeleteSharedArray(eigenocc,'eigenocc')
                 LEND63t = .true.
                 call msg('CalcAk: not all Ak submatrices are done, you can restart from TAPE63t')
                 GOTO 100 
                 end if !prova restart

               do jvirt = 1, nbas(isymvirt)

                  pfcol = 0.0_KREAL
                  virtatom = invbas(jvirt,isymvirt)
                  virtbas = npart(jvirt,isymvirt)

                  do jocc = 1, nbas(isymocc)

                     occatom = invbas(jocc,isymocc)
                     occbas = npart(jocc,isymocc)

                     maxatom = max(virtatom, occatom)
                     minatom = min(virtatom, occatom)

                     ipair = ((maxatom-1)*maxatom)/2 + minatom

                     if ( .not. gPoltddft%skipnb(ipair)) cycle

                     pairlabel = 'atompair'
                     call csaddi (pairlabel, '_', maxatom)
                     call csaddi (pairlabel, '_', minatom)

                     call kfopsc (iu62, pairlabel)

                     basa = atombas(maxatom+1)-atombas(maxatom)
                     basb = atombas(minatom+1)-atombas(minatom)

                     if (maxatom == minatom ) then
                        cplbas = ((basa+1)*basa)/2
                     else
                        cplbas = basa*basb
                     endif

                     allocate (ibasa(cplbas))
                     allocate (ibasb(cplbas))
                     allocate (skip(cplbas+1))
                     ibasa = 0
                     ibasb = 0
                     skip = 0 

                     maxbas = max(virtbas, occbas)
                     minbas = min(virtbas, occbas)

                     call kfread (iu62, 'basa', ibasa)
                     call kfread (iu62, 'basb', ibasb)
                     call kfread (iu62, 'skip', skip)

                     do j = 1,  cplbas
                        if (ibasa(j) ==0 .or. ibasb(j) ==0) EXIT
                        if (maxbas == ibasa(j) .and. minbas == ibasb(j)) then
                    
                           allocate (pairfitkf(skip(j+1)-skip(j)))
                           allocate (idxpf(skip(j+1)-skip(j)))

                           idxpf = 0
                           pairfitkf = 0.0_KREAL

                           call kfopvr(iu62,'idxpf')
                           call kfskvr(iu62, skip(j))
                           call kfrdni(iu62,'%', idxpf,skip(j+1)-skip(j), 1) 

                           call kfopvr(iu62,'pairfit')
                           call kfskvr(iu62, skip(j))
                           call kfrdnr(iu62,'%', pairfitkf, skip(j+1)-skip(j), 1)

                           do l = 1, skip(j+1)-skip(j)
                              pfcol(idxpf(l)) = pfcol(idxpf(l)) + pairfitkf(l)*eigenocc(jocc,i)
                           enddo 

                           deallocate(idxpf, pairfitkf)
                           exit

                        endif
                     enddo
                  
                     deallocate(ibasa, ibasb, skip)

                  enddo
                  call LockSharedArray (pairfit)
                  pairfit(:,jvirt) =  pairfit(:,jvirt) + pfcol(:) 
                  call UnLockSharedArray (pairfit)

               enddo

               call ppcbnr (pairfit, size(pairfit), 'pairfit')

               allocate(Ak(gDims%nsfos,gPoltddft%nvirtmax(i,isymocc,isymvirt)))

               Ak = 0.0_KREAL
               
               call mmulrr(pairfit,gDims%nsfos,nbas(isymvirt), &
                           eigenvirt(1,virt1+1),nbas(isymvirt),gPoltddft%nvirtmax(i,isymocc,isymvirt), &
                           'NO',Ak)

               allocate(Aktmp(gDims%nsfos))
               Aktmp = 0.0_KREAL

               do n = 1, gPoltddft%nvirtmax(i,isymocc,isymvirt) 
                 k = k + 1
                 Aktmp(1:gDims%nsfos) = Ak(1:gDims%nsfos, n)
                 kstr = ' '
                 call csputi (kstr, k)
                 call kfwrnr(iu63t, 'Ak'//trim(kstr), Aktmp, gDims%nsfos, 1) !prova change from t63 to t63t
               end do
          
               counter = counter + gPoltddft%nvirtmax(i,isymocc,isymvirt) !prova restart
               
               if (print_it.eq.1) then
                  write(line,*) counter !prova max number of submatrices
                  line = "Number of submatrices Ak done = " // trim(line) !prova max number of submatrices
                  call msg(trim(line)) !prova max number of submatrices
               endif

               deallocate(Ak, Aktmp)

            enddo
            call DeleteSharedArray(pairfit,'pairfit')
            call DeleteSharedArray(eigenvirt,'eigenvirt')
            exit irrep_

         end do irrep_
      enddo
      call DeleteSharedArray(eigenocc,'eigenocc')
   enddo

   call kfclfl(iu63t) !prova change from t63 to t63t

   deallocate (pfcol, npart, invbas) !prova collected deallocate
   call kfclfl(iu62)
   call kfdlfl('TAPE62')
   call kfclfl(iu21)

   100 CONTINUE

   call timere('pp')
   call timere ('Ak')

write(iuout,*) "start of CalckAk PIER"

End subroutine CalcAk 
