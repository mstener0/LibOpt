Subroutine  Calcdpadqi(nsymr, nrocc, irrepop)

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
   use ADFGlobalInputModule

   implicit none 

   integer(KINT), intent(in)  :: nsymr, irrepop(gSymmetryInfo%nrep)
   integer(KINT), intent(in)  :: nrocc(gExcitSOVars%nsymt,gExcitOpenMod%nds)

   integer(KINT) :: i, j, jj, f, l, k, y, irrep, nrpq, c, d 
   integer(KINT) :: jatyp, ja
   integer(KINT) :: max_nvirtmax, nocc, min_noccmin  
   integer(KINT) :: isym 
   integer(KINT) :: jatom, jjatom, maxatom, minatom
   integer(KINT) :: jbas, jjbas, maxbas, minbas 
   integer(KINT) :: lalo, lahi, maa
   integer(KINT) :: basa, basb, cplbas
   integer(KINT) :: ipair
   integer(KINT) :: iu21, iu62, iu65t 
   integer(KINT) :: atombas(gDims%nnuc+1)
   integer(KINT) :: nbas(nsymr), nbasmax, nmo, fake_noccmin(nsymr,nsymr) 

   character(LCHARS) :: symrep(nsymr)

   integer(KINT), allocatable :: npart(:,:), invbas(:,:)
   integer(KINT), allocatable :: ibasa(:), ibasb(:)
   integer(KINT), allocatable :: skip(:), idxpf(:)

   real(KREAL), allocatable   :: pairfitkf(:)
   real(KREAL), allocatable   :: pfcol(:)

   character(LCHARS)          :: pairlabel

   real(KREAL), allocatable   :: Ak(:,:)

   character(LCHARS)          :: ystr

   real(KREAL), pointer       :: eigen(:,:) 
   real(KREAL), pointer       :: pairfit(:,:)
   character(LCHARS) :: line
   integer(KINT) :: iexist65t, threshold, n_submatrices_restart, kk

   call timers ('Ak')

   call gInput%Get ('PolTDDFT%N_FitOrb', threshold)

   call kfopfl (iu21, gADFFiles%main)
   call kfopvr (iu21, 'Symmetry%nsym')
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

   fake_noccmin = gPoltddft%noccmin
   do c = 1, nsymr
      do d = 1, nsymr
         if (fake_noccmin(c,d) == 0) fake_noccmin(c,d) = 1000000
      enddo
   enddo

   do isym = 1, nsymr
       max_nvirtmax = maxval(gPoltddft%nvirtmax(:,:,isym))
       min_noccmin = minval(fake_noccmin(isym,:))
       nocc = nrocc(isym, 1)
       if (nocc+max_nvirtmax == 0) cycle
   end do
   ipair = 0 

   call timers('pp')

   iexist65t = 0
   if (kfexfl('TAPE65t')) iexist65t = 1
   call ppNodeBarrier
   call ppcbi(iexist65t,'iexst65t')

   n_submatrices_restart = 0

   if (iexist65t /= 0) then
     call kfopflpp(iu65t, 'TAPE65t') 
     call kfclfl(iu65t)

     call kfopfl(iu65t, 'TAPE65t')
     call kfopsc(iu65t, 'MatrixAkt') 
     call kfread(iu65t, 'MatrixAkt%n_submatrices_restart', n_submatrices_restart)
   else
     call kfcrfl(iu65t, 'TAPE65t')
     call kfcrsc(iu65t, 'MatrixAkt')
   end if

   k = n_submatrices_restart
   threshold = threshold + n_submatrices_restart
   allocate (pfcol(gDims%nsfos))
   pfcol = 0.0_KREAL  
   kk = 0
   call kfopfl (iu62, 'TAPE62')

   do isym = 1, nsymr
      y = 0 
      call kfcrsc(iu65t,trim(symrep(isym))//'_occ')
      call kfcrsc(iu65t,trim(symrep(isym))//'_virt')
      max_nvirtmax = maxval(gPoltddft%nvirtmax(:,:,isym))
      min_noccmin = minval(fake_noccmin(isym,:))
      nocc = nrocc(isym, 1)
      if (nocc+max_nvirtmax == 0) cycle
         call kfopsc(iu21, symrep(isym))
         call kfread(iu21, 'nmo_A', nmo)
         call NewSharedArray(eigen,'eigen',(/nbas(isym),nmo/))
         if (AmIOwner(eigen)) call kfread(iu21, 'Eigen-Bas_A', eigen)
         if (IsSharedArray(eigen)) call ppNodeBarrier
         call NewSharedArray(pairfit,'pairfit',(/gDims%nsfos,nbas(isym)/))
         do i = min_noccmin, nocc+max_nvirtmax
               if (AmIOwner(pairfit)) pairfit = 0.0_KREAL
               if (IsSharedArray(pairfit)) call ppNodeBarrier
               if ((kk + 1) <= n_submatrices_restart) then
                  kk = kk +1 
                  cycle
               end if
               if ((k +1)  > threshold) then
                  call kfwrite(iu65t, 'n_submatrices_restart', k)
                  call kfclfl(iu65t)
                  call kfclfl(iu62)
                  call kfdlfl('TAPE62')
                  call kfclfl(iu21)
                  if (IsSharedArray(pairfit)) call ppNodeBarrier
                  call stopit ('Blocco da finire Ak')
               end if

               do j = 1, nbas(isym)
                  pfcol = 0.0_KREAL
                  jatom = invbas(j,isym)
                  jbas = npart(j,isym)

                  do jj = 1, nbas(isym)

                     jjatom = invbas(jj,isym)
                     jjbas = npart(jj,isym)

                     maxatom = max(jatom, jjatom)
                     minatom = min(jatom, jjatom)

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

                     maxbas = max(jbas, jjbas)
                     minbas = min(jbas, jjbas)

                     call kfread (iu62, 'basa', ibasa)
                     call kfread (iu62, 'basb', ibasb)
                     call kfread (iu62, 'skip', skip)
                        do f = 1,  cplbas
                           if (ibasa(f) ==0 .or. ibasb(f) ==0) EXIT
                           if (maxbas == ibasa(f) .and. minbas == ibasb(f)) then
                              allocate (pairfitkf(skip(f+1)-skip(f)))
                              allocate (idxpf(skip(f+1)-skip(f)))

                              idxpf = 0
                              pairfitkf = 0.0_KREAL
                              call kfopvr(iu62,'idxpf')
                              call kfskvr(iu62, skip(f))
                              call kfrdni(iu62,'%', idxpf,skip(f+1)-skip(f), 1) 

                              call kfopvr(iu62,'pairfit')
                              call kfskvr(iu62, skip(f))
                              call kfrdnr(iu62,'%', pairfitkf, skip(f+1)-skip(f), 1)
                                 do l = 1, skip(f+1)-skip(f)
                                    if (i == 0) write(iuout,*) 'i = 0'
                                    pfcol(idxpf(l)) = pfcol(idxpf(l)) + pairfitkf(l)*eigen(jj,i)
                                 enddo 
                              deallocate(idxpf, pairfitkf)
                              exit

                           endif
                        enddo
                     deallocate(ibasa, ibasb, skip)
                  enddo
                  call LockSharedArray (pairfit)
                  pairfit(:,j) =  pairfit(:,j) + pfcol(:) 
                  call UnLockSharedArray (pairfit)
               enddo
         call ppcbnr (pairfit, size(pairfit), 'pairfit')

         allocate(Ak(gDims%nsfos,1)) 

         Ak = 0.0_KREAL


         call mmulrr(pairfit,gDims%nsfos,nbas(isym), &
                     eigen(1,i),nbas(isym),1,'NO',Ak)
         k = k + 1
         if (i > nocc) then
            y = i - nocc
            ystr = ' '
            call csputi (ystr, y)
            call kfwrnr(iu65t,trim(symrep(isym))//'_virt%dens_fit_'//trim(ystr),Ak(1:gDims%nsfos,1),gDims%nsfos,1)
         else
             y = i - min_noccmin + 1
            ystr = ' '
            call csputi (ystr, y)
            call kfwrnr(iu65t,trim(symrep(isym))//'_occ%dens_fit_'//trim(ystr),Ak(1:gDims%nsfos,1),gDims%nsfos,1)
         end if

         deallocate(Ak)

      enddo
      call DeleteSharedArray(pairfit,'pairfit')
      call DeleteSharedArray(eigen,'eigen')
   enddo

   call timere('pp')

   call kfclfl(iu65t) !prova cambiato da t63 a t63t

   deallocate (pfcol, npart, invbas) !prova riunito deallocate

   call kfclfl(iu62)
   call kfclfl(iu21)

   call timere ('Ak')

write(iuout,*) "end of Calcdpq... PIER"

End subroutine Calcdpadqi
