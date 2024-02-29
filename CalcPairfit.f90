subroutine CalcPairfit ()

   use KF
   use DependencyInfo
   use DimensionsModule
   use SpinsuffixModule
   use ConstModule
   use FitsymModule
   use ScfdataModule
   use MasterModule
   use ModelDataModule
   use FactorialsModule
   use PoltddftModule
   use ADFFilesModule

   use ppCommonsModule
   use ADFGlobalInputModule
   use SharedArraysModule
   use SharedArraysUtilModule

   use Vartypes
   implicit real(KREAL) (a-h, o-z)
   implicit integer(KINT) (i-n)

   integer(KINT) :: iskf(4,gDims%niskf), na1ptr(gDims%nnuc+1)
   integer(KINT), allocatable :: numcom(:)
   real(KREAL), allocatable   :: cofcom(:), tfit(:)  

   real(KREAL), allocatable   :: sfit(:), sfitab(:,:)

   real(KREAL), allocatable   :: pairfit(:,:), pairfit3c(:,:)
   real(KREAL), allocatable   :: tfita(:,:), tfitb(:,:)

   integer(KINT), allocatable :: basa(:), basb(:)
   integer(KINT), allocatable :: skip(:)
   integer(KINT), allocatable :: ibasa(:), ibasb(:)

   real(KREAL), allocatable   :: descale(:), scale(:)
   integer(KINT), allocatable :: ind(:), idxpf(:)

   real(KREAL), allocatable :: cfitints(:), amatfit(:,:), vectfit(:)
   integer, allocatable :: indea(:), indeb(:)

   integer(KINT) :: inom(496), minp(31), itpow(31)

   integer(KINT), parameter   :: ncxamax = 2892
   real(KREAL)                :: cxa(ncxamax), vec(3)
   real(KREAL), pointer       :: intfit(:)
   real(KREAL), pointer       :: sfitfull(:,:)
   real(KREAL), allocatable   :: packed(:)
   integer(KINT)              :: nx

!  ====================================================================
!  purpose: write pair fitting integrals
!
!  ====================================================================

   logical :: a1fit, centr1, lofi, altRhoF

   parameter ( one = 1d0,       &
                 int1 = 1, int0 = 0)
   real(KREAL), parameter     :: thres = 1.0e-12_KREAL

   integer(KINT) :: i, l
   integer(KINT) :: nfcom_dummy

!  ---------------
!  dynamic storage
!  ---------------

   integer(KINT)  :: jasym(gDims%npeq)
   integer(KINT)  :: ja1ok(gDims%npeq)
   integer(KINT)  :: neighsym(gDims%nnuc*gDims%nnuc)
   integer(KINT)  :: nprs2nts(gDims%nnuc*gDims%nnuc)
   integer(KINT)  :: nnbr(gDims%nnuc)

   character(LCHARS)          :: pairlabel

!  ===============================
!  initialization / print switches
!  ===============================

   allocate (gPoltddft%skipnb(gDims%nsfos*(gDims%nsfos+1)/2))
   gPoltddft%skipnb = .false.


   allocate (numcom(1:gDims%na1cof))
   numcom = 0
   allocate (cofcom (gDims%na1cof))
   cofcom = 0.0_KREAL

   ltfit = max(gDims%lrl,gDims%maxsf)

   call timers('Pairfit')

   call timers('read')

!  ------------------
!  read symmetry info
!  ------------------

   call gInput%Get ('altRhoF', altRhoF)
   
   call kfopfl (iu, gADFFiles%main)
   if (gFitsym%lnosft.and..not.altRhoF) then
      call kfopsc (iu, 'SymFit')
   else
      call kfopsc (iu, 'Symmetry')
   end if
   call kfread (iu, 'jasym', jasym)
   call kfread (iu, 'ja1ok', ja1ok)
   call kfclfl (iu)

!  -----------------------------
!  open tape with fit integrals:
!  iuf: fit integrals themselves
!  -----------------------------

   call kfopfl (iuf, 'TAPE11')

   call rdfa1c (.false., gDims%nnuc, gDims%nsetat, gDims%niskf, gDims%na1cof, gNuclei%nratst, &
                gNuclei%noat, nfcom_dummy, na1ptr, iskf, numcom, cofcom)
  
!  ----------------------------
!  Read neighbor info for atoms
!  ----------------------------

   call kfopfl (iunbr, 'TAPE15')
   call kfopsc (iunbr, 'Symmetry Neighbours')
   call kfopvr (iunbr, 'Num of neighs')
   call kfopfl (iunei, 'TAPE15')
   call kfopsc (iunei, 'Symmetry Neighbours')
   call kfopvr (iunei, 'Neighbours')
   call kfopfl (iunts, 'TAPE15')
   call kfopsc (iunts, 'Symmetry Neighbours')
   call kfopvr (iunts, 'Sym pair number')
!
   call kfrdni (iunbr, '%', nnbr, gDims%nnuc, 1)

!  ------------------------------------
!  calculate total number of neighbours
!  ------------------------------------

!write(*,*) 'nnbr', nnbr

   ntotbr = 0
   do inuc = 1, gDims%nnuc
      ntotbr = ntotbr + nnbr(inuc)
   end do
!
   call kfrdni (iunei, '%', neighsym, ntotbr, 1)
   call kfrdni (iunts, '%', nprs2nts, ntotbr, 1)
!
   call kfclfl (iunbr)
   call kfclfl (iunei)
   call kfclfl (iunts)


   call kfopfl (iud, 'TAPE11')
   call kfopsc (iud, 'Fitint')
   call kfread (iud, 'Fitint%clen', lirec)
   allocate (cfitints(lirec),stat=istat)
   cfitints =0.0_KREAL
   call chckmem (0,'fitint.d90:cfitints',istat)
   call kfrdnr (iud, 'Fitint%cfitints', cfitints, lirec, 1) 

   call kfclsc(iud)   

   call kfclfl (iud)

   call timere('read')

!  =====================
!  start loop over spins
!  =====================

!  --------------------------
!  position fit data variable
!  --------------------------

   call kfopfl(ius,'TAPE11')
   call kfopnp(ius, 'Fitdata%sfac')


   call kfopnp (iuf, 'FitIntegrals%data')
!
!   call sicsfitreset (.true., .false., lsicsfiton)

!  -----------------------------------------------
!  read total sfit to calculate 3 center integrals
!  -----------------------------------------------
!

   call kfopfl (iusf, 'TAPE58')

   nx = gDims%nsfos*(gDims%nsfos+1)/2  

   call NewSharedArray(sfitfull,'matrix%sfitfull',(/gDims%nsfos,gDims%nsfos/))
   call SetSharedArrayToZero(sfitfull)
   if (AmIOwner(sfitfull)) then
      allocate(packed(nx))
      packed = 0.0_KREAL
      call kfrdnp (iusf, 'Fitsfit%sfit', packed, nx, 1)
      call mcoppl(gDims%nsfos, packed, sfitfull)
      call mcoplu (gDims%nsfos, sfitfull)
      deallocate(packed)
   end if
   call ppNodeBarrier

   call kfclfl(iusf)

!  -------------------------------
!  initialize inom, minp and itpow 
!  -------------------------------

   call nom (inom, minp, itpow)

   lfna = 0
   do ityp = 1, gDims%ntyp
      lfna = max(lfna,2*gTypesc%nftes(ityp))
   end do

   allocate (ind(lfna))
   ind =0

   allocate(intfit(gDims%nsfos))
   intfit = 0.0_KREAL   

   call kfcrfl (iup, 'TAPE62')
   call kfclfl (iup)

   kkx  = 0
   nts  = 1

   call timers ('pp')

!  -----------------------------------------------
!  s t a r t   l o o p   o v e r   a t o m s   (a)
!  -----------------------------------------------

   ibrdone = 0
   jatyp_: do jatyp = 1, gDims%ntyp
!
      lalo = gTypesc%nbptr(jatyp)
      lahi = gTypesc%nbptr(jatyp+1) - 1
      maa  = lahi - lalo + 1
      nalo = gTypesc%nbaspt(jatyp)
      nahi = gTypesc%nbaspt(jatyp+1) - 1

      lmxatyp = gFittyp%lqfitx(jatyp) + gBastyp%lqbasx(jatyp)

      do ja = gTypesc%nqptr(jatyp), gTypesc%nqptr(jatyp+1) - 1

         kalo = na1ptr(ja)
         kahi = na1ptr(ja+1) - 1

         do la = kalo, kahi
            ka = iskf(1,la)

            lxa       = gFittyp%lqfit(ka)
            tlpa      = lxa + lxa + 1
            nra       = gFittyp%nqfit(ka) + 2
            alfa      = gFittyp%alffit(ka)

            if (lxa == 0) then
               sa = gFactorials%factqu(nra-1,1)/alfa**(nra)
               intfit(la) = gConstants%twopi*(sa+sa)/tlpa
            end if
         enddo

         if (maa > 0) then
             kkx = gNuclei%indxbs(ja)
         end if

!        --------------------------------------------------------------
!        s t a r t   i n n e r  l o o p   o v e r   n e i g h b o u r s
!        --------------------------------------------------------------

         nsymbr = nnbr(ja)

         inbr_: do inbr = 1, nsymbr

!           ----------------------------------------------
!           atom number and number of symmetry unique pair
!           ----------------------------------------------

            jb    = neighsym(ibrdone+inbr) 
    
            nts   = nprs2nts(ibrdone+inbr)
            jbtyp = gNuclei%iat2ityp(jb)

            lblo = gTypesc%nbptr(jbtyp)
            lbhi = gTypesc%nbptr(jbtyp+1) - 1
            mbb  = lbhi - lblo + 1
            nblo = gTypesc%nbaspt(jbtyp)
            nbhi = gTypesc%nbaspt(jbtyp+1) - 1
            if (maa==0 .or. mbb==0) cycle inbr_
!
            llx = gNuclei%indxbs(jb)

            kblo = na1ptr(jb)
            kbhi = na1ptr(jb+1) - 1

!           -----------------------------
!           pair number and pair distance
!           -----------------------------

            dx   = gNuclei%xyznuc(1,ja) - gNuclei%xyznuc(1,jb)
            dy   = gNuclei%xyznuc(2,ja) - gNuclei%xyznuc(2,jb)
            dz   = gNuclei%xyznuc(3,ja) - gNuclei%xyznuc(3,jb)
            rab2 = dx**2 + dy**2 + dz**2
            rab  = sqrt(rab2)
!
            centr1 = ja==jb
!
            a1fit = ja1ok(nts)==1
            if (a1fit) then
               nfaa = na1ptr(ja+1) - na1ptr(ja)
               nfbb = na1ptr(jb+1) - na1ptr(jb)
            else
               nfaa = gTypesc%nftes(jatyp)
               nfbb = gTypesc%nftes(jbtyp)
            end if
            if (centr1) then
               nfab = nfaa
            else
               nfab = nfaa + nfbb
            end if

            ipair = (ja-1)*ja/2+jb
            gPoltddft%skipnb(ipair) = .true.
!
!           ---------------
!           one center case
!           ---------------

            if (centr1) then

!              ----------------------------------------
!              determine normalization of fit integrals
!              ----------------------------------------

               allocate (sfit(gDims%maxsf))

               sfit = 0.0_KREAL

               call fa1soc (lfna, int0, kalo, kahi, iskf, ind, sfit)

               allocate (descale(nfab))
               descale = 0.0_KREAL             
 
               idx = 0
               do ifit = 1, nfab
                  idx         = idx + ifit
                  descale(ifit) = sqrt(sfit(idx))
               end do

               mxmx = gBasfcn%ndex(nfab+1)
               call kfrdnp(ius, '%', sfit, mxmx, int1)

               deallocate (sfit)

!              -------------------------
!              sfit(i) = rho(aa) * fa(i)
!              -------------------------

               kk   = kkx
               ipmt = 0

               allocate(basa(maa*(maa+1)/2))
               allocate(basb(maa*(maa+1)/2))
               basa = 0
               basb = 0

               allocate(pairfit(gDims%nsfos,maa*(maa+1)/2))
               pairfit = 0.0_KREAL    

               allocate(tfit(nfab))
               tfit = 0.0_KREAL     

               call kfopfl (iup, 'TAPE62')

               pairlabel = 'atompair'
               call csaddi (pairlabel, '_', ja)
               call csaddi (pairlabel, '_', ja)

               call kfcrsc (iup, pairlabel)

               call kfcrni(iup,'basa',maa*(maa+1)/2)
               call kfcrni(iup,'basb',maa*(maa+1)/2)
               call kfcrni(iup,'skip',maa*(maa+1)/2+1)

               i = 1

               do la = lalo, lahi
                  kk  = kk + 1
                  ndx = gBasfcn%ndex(kk) + llx
                  ll  = llx
                  do lb = lblo, la
                     ll   = ll + 1
                     ndx  = ndx + 1
                     ipmt = ipmt + 1

                     call kfrdnp (iuf, '%', tfit, nfab, int1)

                     basa(i) = gNuclei%indxbs(ja) + la - lalo + 1
                     basb(i) = gNuclei%indxbs(jb) + lb - lblo + 1

                     pairfit(kalo:kahi,i) = tfit(1:nfab) * descale(1:nfab)
               
                     i = i + 1

                  end do
               end do

               deallocate(descale)
               deallocate(tfit)

!              ----------------
!              cycle for fbbaba
!              ----------------

               do jb = 1, gDims%nnuc 
 
                  if (jb == ja) cycle
                  
                  jbtyp = gNuclei%iat2ityp(jb)
 
                  lblo = gTypesc%nbptr(jbtyp)
                  lbhi = gTypesc%nbptr(jbtyp+1) - 1
                  mcc  = lbhi - lblo + 1
                  if (maa==0 .or. mcc==0) cycle 
!
                  kblo  = na1ptr(jb)
                  kbhi  = na1ptr(jb+1) - 1

                  lmxbtyp = gFittyp%lqfitx(jbtyp) + gBastyp%lqbasx(jbtyp)
                  lmxabbas = gBastyp%lqbasx(jatyp) + gBastyp%lqbasx(jbtyp)

!                 -----------------------------
!                 pair number and pair distance
!                 -----------------------------

                  vec(1)   = gNuclei%xyznuc(1,jb) - gNuclei%xyznuc(1,ja)
                  vec(2)   = gNuclei%xyznuc(2,jb) - gNuclei%xyznuc(2,ja)
                  vec(3)   = gNuclei%xyznuc(3,jb) - gNuclei%xyznuc(3,ja)

                  lmxab = max(lmxatyp,lmxbtyp,lmxabbas)
                  call mkxyzrotmat (vec, lmxab, ncxamax, ncxa, cxa, rab)
 
                  a1fit = ja1ok(nts)==1
                  if (a1fit) then
                     nfbb = na1ptr(jb+1) - na1ptr(jb)
                  else
                     nfbb = gTypesc%nftes(jbtyp)
                  end if

                  irec = nfbb*(maa*(maa+1)/2)

                  allocate (tfit(irec))
                  tfit = 0.0_KREAL

                  call fbbaba(nfbb, lalo, lahi, kblo, kbhi, iskf, numcom, cofcom,   &
                              rab, inom, minp, itpow, .true., ncxa, cxa, nalo,      &
                              irec, tfit)

                  i = 1
                  j = 1  

                  do la = lalo, lahi
                     do lb = lalo, la
                        pairfit(kblo:kbhi,i) = tfit(j:j+nfbb-1)  

                        j = j + nfbb
                        i = i + 1
                     enddo
                  enddo

                  deallocate (tfit)

               enddo  

               allocate(ibasa(maa*(maa+1)/2))
               allocate(ibasb(maa*(maa+1)/2))
               ibasa = 0
               ibasb = 0

               allocate (skip ((maa*(maa+1)/2)+1))
               skip = 0

               l = 0
               do i = 1, maa*(maa+1)/2
                  if (any(ABS(pairfit(1:gDims%nsfos,i)) > thres )) then
                     l = l + 1 
                     skip(l+1) =skip(l)
                     do j = 1, gDims%nsfos
                        if (ABS(pairfit(j,i))>thres) skip(l+1) = skip(l+1) + 1
                     enddo
                     ibasa(l) = basa(i)
                     ibasb(l) = basb(i)

                     icpl = skip(l+1)
                  endif
               enddo

               call kfopvr (iup, 'basa')
               call kfwrni (iup, '%', ibasa, maa*(maa+1)/2, 1)
               call kfopvr (iup, 'basb')
               call kfwrni (iup, '%', ibasb, maa*(maa+1)/2, 1)
               call kfopvr (iup, 'skip')
               call kfwrni (iup, '%', skip, maa*(maa+1)/2+1, 1)

               deallocate (skip)

               deallocate (basa, basb)

               deallocate (ibasa, ibasb)
       
               allocate (idxpf(icpl))
               idxpf = 0

               call kfcrni(iup,'idxpf',icpl)
               call kfcrnr(iup,'pairfit',icpl) 

               l = 0
               do i = 1, maa*(maa+1)/2
                  if (any(ABS(pairfit(1:gDims%nsfos,i)) > thres )) then
                     do j = 1, gDims%nsfos
                        if (ABS(pairfit(j,i))>thres) then
                           l = l + 1 
                           idxpf(l) = j
                           call kfwrnr (iup, '%', pairfit(j,i), 1, 1) 
                        endif
                     enddo
                  endif
               enddo

               call kfopvr (iup, 'idxpf')
               call kfwrni (iup, '%', idxpf, icpl, 1)
         
               call kfclfl (iup)

               deallocate (idxpf)

               deallocate (pairfit)

!           ---------------
!           two center case
!           ---------------

            else

               allocate (sfit(gDims%maxsf))
               sfit = 0.0_KREAL

               lmxbtyp = gFittyp%lqfitx(jbtyp) + gBastyp%lqbasx(jbtyp)
               lmxabbas = gBastyp%lqbasx(jatyp) + gBastyp%lqbasx(jbtyp)

               vec(1)   = gNuclei%xyznuc(1,jb) - gNuclei%xyznuc(1,ja)
               vec(2)   = gNuclei%xyznuc(2,jb) - gNuclei%xyznuc(2,ja)
               vec(3)   = gNuclei%xyznuc(3,jb) - gNuclei%xyznuc(3,ja)

               lmxab = max(lmxatyp,lmxbtyp,lmxabbas)
               call mkxyzrotmat (vec, lmxab, ncxamax, ncxa, cxa, rab)

               call fa1soc (lfna, int0, kalo, kahi, iskf, ind, sfit)
               call fa1soc (lfna, nfaa, kblo, kbhi, iskf, ind, sfit)
               call fa1stc (nfaa, kalo, kahi, kblo, kbhi, iskf, sfit, numcom, cofcom, ncxa, cxa, &
                            rab, inom, minp, itpow)

               allocate (descale(nfab))
               allocate (scale(nfab))
               descale = 0.0_KREAL
               scale = 0.0_KREAL

               idx = 0
               do ifit = 1, nfab
                  idx         = idx + ifit
                  descale(ifit) = sqrt(sfit(idx))
                  scale(ifit) = one/sqrt(sfit(idx))
               end do

!-----------------------------------------
! 3 center 30/05/17
!-----------------------------------------

               nrig = nfaa+nfbb
               ncolon = maa*mbb
               allocate (amatfit(nrig,ncolon), vectfit(nrig),indea(ncolon),indeb(ncolon) )
               indea = 0
               indeb = 0
               amatfit = 0.0_KREAL
               vectfit = 0.0_KREAL
               lofi = .true.

               i0 = 0

               call fa1ttc3c (nfaa, lalo, lahi, lblo, lbhi, kalo, kahi, iskf, numcom, cofcom,   &
                             rab, inom, minp, itpow, .true., scale,  &
                             ncxa, cxa, nalo, nblo, amatfit, nrig, ncolon, i0, indea, indeb)

!              ---------------------------
!              fit integral f(b) b(B) b(A)
!              ---------------------------
               
               i0 = nfaa
               call fa1ttc3c (nfbb, lblo, lbhi, lalo, lahi, kblo, kbhi, iskf, numcom, cofcom,   &
                             rab, inom, minp, itpow, .false., scale(nfaa+1),          &
                             ncxa, cxa, nblo, nalo, amatfit, nrig, ncolon, i0, indea, indeb )

               sfit =0.0_KREAL
               mxmx = gBasfcn%ndex(nfab+1)
               call kfrdnp(ius, '%', sfit, mxmx, int1)

               if (lofi) then
                  if (nfab.ne.nrig) stop ' nrig '
                  do icol = 1, ncolon
                     vectfit = 0.0_KREAL
                     call msolpp (nfab, sfit, amatfit(1,icol), vectfit)
                     amatfit(1:nfab,icol) = vectfit(1:nfab)
                  enddo

                  do irig = 1, nfab
                     amatfit(irig,:) = amatfit(irig,:)*scale(irig)
                  enddo

                  lofi = .false.
               endif

               deallocate (scale)
               deallocate (sfit)
               deallocate(vectfit, indea, indeb)

               allocate (sfitab(gDims%nsfos,nfab))
               sfitab = 0.0_KREAL

               sfitab(1:gDims%nsfos,1:nfaa) = sfitfull(1:gDims%nsfos,kalo:kahi)
               sfitab(1:gDims%nsfos,nfaa+1:nfab) = sfitfull(1:gDims%nsfos,kblo:kbhi)
               sfitab(kalo:kahi, 1:nfab) = 0.0_KREAL
               sfitab(kblo:kbhi, 1:nfab) = 0.0_KREAL

               allocate (pairfit3c(gDims%nsfos,ncolon))
               pairfit3c = 0.0_KREAL
 
               call mmulrr(sfitab, gDims%nsfos, nfab, amatfit, nfab, ncolon, 'NO', pairfit3c)
                        
               deallocate(amatfit, sfitab)

!-------------------------------------------

               call kfopfl (iup, 'TAPE62')

               pairlabel = 'atompair'
               call csaddi (pairlabel, '_', ja)
               call csaddi (pairlabel, '_', jb)

               call kfcrsc (iup, pairlabel)
         
               allocate(basa(maa*mbb))
               allocate(basb(maa*mbb))
               basa = 0
               basb = 0

               allocate (tfita(nfaa,maa*mbb))
               tfita = 0.0_KREAL

               icpl = 0
               ina = 0

               do na = nalo, nahi
                  nfa    = gBasfcn%ndex(gBastyp%lqbas(na)+2)
           
                  inb = 0
                  do nb = nblo, nbhi
                     nfb = gBasfcn%ndex(gBastyp%lqbas(nb)+2)

!                    -------------------------------------------
!                    skip if basis functions are non-overlapping
!                    -------------------------------------------

                     if ((gBastyp%bradfit(na)+gBastyp%bradfit(nb))<rab) then 
                        do j = 1, nfb
                           do m = 1, nfa
                              icpl = icpl + 1

                              basa(icpl) = gNuclei%indxbs(ja) + ina + m
                              basb(icpl) = gNuclei%indxbs(jb) + inb + j
                           enddo
                        enddo  
                        inb = inb + nfb
                        cycle
                     endif

                     irec = nfaa*nfa*nfb

                     allocate(tfit(irec))
                     tfit = 0.0_KREAL

                     call kfrdnp (iuf, '%', tfit, irec, int1)

                     l = 1

                     do j = 1, nfb
                        do m = 1, nfa
!
                           icpl = icpl + 1

                           basa(icpl) = gNuclei%indxbs(ja) + ina + m
                           basb(icpl) = gNuclei%indxbs(jb) + inb + j

                           tfita(1:nfaa,icpl) = tfit(l:l+nfaa-1) * descale(1:nfaa)
                           l = l + nfaa
                        enddo
                     enddo

                     deallocate(tfit)
                     inb = inb + nfb
                  end do 
                  ina = ina + nfa
               end do

               call kfcrni(iup, 'basa', maa*mbb)
               call kfcrni(iup, 'basb', maa*mbb)
               call kfcrni (iup, 'skip', maa*mbb+1)

               allocate (tfitb(nfbb,maa*mbb))
               tfitb = 0.0_KREAL

               i = 1

               do nb = nblo, nbhi
                  nfb = gBasfcn%ndex(gBastyp%lqbas(nb)+2)
!        
                  l1 = i

                  nfap = 1

                  do na = nalo, nahi
                     nfa = gBasfcn%ndex(gBastyp%lqbas(na)+2)

!                    -------------------------------------------
!                    skip if basis functions are non-overlapping
!                    -------------------------------------------

                     if ((gBastyp%bradfit(na)+gBastyp%bradfit(nb))<rab) then

                        if(nfa/=nfap ) then
                          l1 = l1 + ((nfa-nfap)*(i-1))
                          nfap = nfa
                        endif
                        do j = 1, nfa
                           if (j > 1) l1 = l1 + 1 - (nfb-1)*nfa
                           do m = 1, nfb
                              if (m > 1) l1 = l1 + nfa
                           enddo
                        enddo
                        l1 = l1 + (mbb*nfa) -nfa +1
                        l1 = l1 - (nfb-1)*nfa
                        cycle
                     endif

                     irec = nfa*nfb*nfbb

                     allocate(tfit(irec))
                     tfit = 0.0_KREAL

                     call kfrdnp (iuf, '%', tfit, irec, 1)

                     l2 = 1

                     if (nfa/=nfap ) then
                       l1 = l1 + ((nfa-nfap)*(i-1))
                       nfap =nfa
                     endif

                     do j = 1, nfa
                        if (j > 1) l1 = l1 + 1 - (nfb-1)*nfa 
                        do m = 1, nfb
                           if (m > 1) l1 = l1 + nfa 

                           tfitb(1:nfbb,l1) = tfit(l2:l2+nfbb-1) * descale(nfaa+1: nfab)       

                           l2 = l2 + nfbb
                        enddo                        
                     enddo

                     deallocate(tfit)
                     l1 = l1 + (mbb*nfa) -nfa +1 
                     l1 = l1 - (nfb-1)*nfa

                  end do 
                  i = i + nfb

               end do 

               deallocate (descale)
               allocate (pairfit(gDims%nsfos,1))
               pairfit = 0.0_KREAL

               allocate(ibasa(maa*mbb))
               allocate(ibasb(maa*mbb))
               ibasa = 0
               ibasb = 0

               allocate (skip(maa*mbb+1))
               skip = 0

               icpl = 0
               l = 0

               do i = 1,maa*mbb
                  if (basa(i) == 0 .or. basb(i) ==0) cycle
                  
                  pairfit(1:gDims%nsfos,1) = pairfit3c(1:gDims%nsfos,i)
                  pairfit(kblo:kbhi,1) = tfitb (1:nfbb,i)
                  pairfit(kalo:kahi,1) = tfita (1:nfaa,i)

                  if (any(ABS(pairfit(1:gDims%nsfos,1)) > thres )) then
                     l= l+1
                     skip(l+1) =skip(l)
                     do j = 1, gDims%nsfos
                        if (ABS(pairfit(j,1))>thres) skip(l+1) = skip(l+1) + 1
                     enddo
                     ibasa(l) = basa(i)
                     ibasb(l) = basb(i)

                     icpl = skip(l+1)
                  endif

               enddo

               call kfopvr (iup, 'basa')
               call kfwrni (iup, '%', ibasa, maa*mbb, 1)
               call kfopvr (iup, 'basb')
               call kfwrni (iup, '%', ibasb, maa*mbb, 1)
               call kfopvr (iup, 'skip')
               call kfwrni (iup, '%', skip(1:maa*mbb+1), maa*mbb+1, 1)

               deallocate (skip)
               deallocate (ibasa, ibasb)
               allocate (idxpf(icpl))
               idxpf =0

               call kfcrni (iup, 'idxpf', icpl) 
               call kfcrnr (iup, 'pairfit', icpl)
               call kfopvr (iup, 'pairfit')

               l = 0

               do i = 1,maa*mbb
                  if (basa(i) == 0 .or. basb(i) ==0) cycle

                  pairfit(1:gDims%nsfos,1) = pairfit3c(1:gDims%nsfos,i)
                  pairfit(kblo:kbhi,1) = tfitb (1:nfbb,i)
                  pairfit(kalo:kahi,1) = tfita (1:nfaa,i)

                  if (any(ABS(pairfit(1:gDims%nsfos,1)) > thres )) then
                     do j = 1, gDims%nsfos
                        if (ABS(pairfit(j,1))>thres) then
                           l = l + 1
                           idxpf(l) = j         
                           call kfwrnr(iup, '%', pairfit(j,1), 1, 1)
                        endif
                     enddo
                  endif

               enddo

               call kfopvr(iup, 'idxpf')  
               call kfwrni(iup, '%', idxpf, icpl, 1)

               call kfclfl(iup)

               deallocate (pairfit3c)
               deallocate (pairfit)
               deallocate (basa, basb)
               deallocate (tfita, tfitb)
               deallocate (idxpf)

            end if                   

!        ------------------------
!        end of the parallel part
!        ------------------------

         end do inbr_

!        ------------------------------------------
!        e n d  l o o p   o v e r   a t o m s   (b)
!        ------------------------------------------

         ibrdone = ibrdone + nsymbr

      end do
   end do jatyp_

   call DeleteSharedArray(sfitfull,'matrix%sfitfull')

   call timere ('pp')

   deallocate (ind)
   deallocate (cofcom)
   deallocate (numcom)

   deallocate (cfitints,stat=istat)
   call chckmem (1,'fitint.d90:cfitints',istat)

!  -----------
!  close files
!  -----------

   call kfopfl (iup, 'TAPE11')
   call kfcrsc (iup, 'intfit')
   call kfwrnr (iup, 'intfit', intfit, gDims%nsfos, 1)

   call kfclfl (iup)

   deallocate(intfit)

   call kfclfl (ius)
   call kfclfl (iuf)

   call timere('Pairfit')

write(iuout,*) "end of CalcPairfit PIER"

end subroutine CalcPairfit
