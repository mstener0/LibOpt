subroutine fa1ttc3c (nfaa, lalo, lahi, lblo, lbhi, kalo, kahi, iskf, numcom, cofcom, rab,       &
                   inom, minp, itpow, fopa, finrm, ncxa, cxa,    &
                   nalo, nblo, amatfit, nrig, ncolon, i0, indea, indeb)

   use KF
   use ConstModule
   use FactorialsModule
   use MasterModule
   use DimensionsModule
   use BitopsInterfaces
!   use FitIntBuffer

   use Vartypes
   implicit real(KREAL) (a-h, o-z)
   implicit integer(KINT) (i-n)

!*if cio
!*copy cexternals
!*ifend

   integer(KINT), intent(in)  :: inom(496), minp(31), itpow(31), iskf(4,gDims%niskf),              &
          numcom(gDims%na1cof)
   logical, intent(in)        :: fopa
   real(KREAL), intent(in)    :: cxa(ncxa), cofcom(gDims%na1cof), finrm(*)
   real(KREAL), intent(inout) :: amatfit(nrig, *)
   integer                    :: indea(ncolon), indeb(ncolon), perm (ncolon)

!  ======================================================================
!  purpose: calculate the two center fit integrals using
!           A1 fit functions
!
!  input: 
!         nfaa - number of fitfunctions on center A
!         lalo - lower index for basis functions on center A
!         lahi - upper index for basis functions on center A
!         lblo - lower index for basis functions on center B
!         lbhi - upper index for basis functions on center B
!         kalo - lower index for A1 fit functions
!         kblo - upper index for A1 fit functions
!         iskf - iskf - pointer array for atomic parts of fit A1-comb's
!         numcom, cofcom - specification of the A1 combination
!         cxa - rotation matrix for two center integrals
!         rab - distance from center A to center B
!         inom - pascal's triangle
!         minp - powers of (-1)
!         itpow - powers of 2
!         fopa - (logical) fit function on center A
!         finrm - normalization factors for the fit functions
!         nalo - lower index in sets of basis functions on a
!         nblo - lower index in sets of basis functions on b
!
!  remark: the calculated fit integrals are written to the current
!          variable of unit as packed reals, thus possibly extending
!          that variable
!
!  called from: fitint
!  ======================================================================

   real(KREAL), parameter   :: zero = 0.0_KREAL
   integer(KINT), parameter :: nctemp = 100, nsinth = 360, nklist = 120, lklist = 7

   integer(KINT) :: klist(3,nklist)
   logical       :: fopb, firstk
   real(KREAL)   :: ctemp(nctemp), sinth(nsinth), sint(nsinth), tfit(gDims%lrl)

   save firstk,klist

   data firstk /.true./

!*if manytimers
!   call timers ('fa1ttc')
!*ifend

#ifdef D_SCM_MANYTIMERS
   call timers ('fa1ttc')
#endif

   if (firstk) then
      firstk = .false.
      call mkxyzlist (lklist, nklist, k, klist)
      if (k/=nklist) call stopit ('klist error. FA1TTC')
   end if

   m0 = 0

   fopb   = .not.fopa
   incidx = 1

!  ------------------------------------------
!  loop over the sets of basis functions on A
!  ------------------------------------------

   k1 = 0
   iseta = nalo - 1
   la    = lalo
   100 if (la>lahi) goto 220
   iseta = iseta + 1

   lxyza = gBastyp%kx(la) + gBastyp%ky(la) + gBastyp%kz(la)
   nfa   = gBasfcn%ndex(lxyza+2)

!  ------------------------------------------
!  loop over the sets of basis functions on B
!  ------------------------------------------

   isetb = nblo - 1
   lb    = lblo
   110 if (lb>lbhi) goto 210
   isetb = isetb + 1

   lxyzb = gBastyp%kx(lb) + gBastyp%ky(lb) + gBastyp%kz(lb)
   nfb   = gBasfcn%ndex(lxyzb+2)
   lrb   = gBastyp%kr(lb)
   alfb  = gBastyp%alf(lb)
   call mkxyzpointers (lxyzb, 0, 0, k0, k, icxb)

!  ------------------------------------
!  make sure the buffer is large enough
!  ------------------------------------

   if (nfaa*nfa*nfb>gDims%lrl) call stopit ('FA1TTC3C: lrl too small')

   if (fopa) incidx = nfb

   llr   = -1
   nprod = nfa*nfb
   if (nprod>nctemp) call stopit ('overflow CTEMP. FA1TTC')

!  -----------------------------------
!  loop over the A1 fit functions on A
!  -----------------------------------

   iast = 0
   ia_: do ia = kalo, kahi
      iast = iast + 1

      ka     = iskf(1,ia)
      kpt    = gCartes%loffst(gFittyp%lqfit(ka))
      ida    = iskf(2,ia)
      naterm = iskf(3,ia)

!     --------------------------------------------------------
!     combine fit function with basis function characteristics
!     --------------------------------------------------------

      lt   = lxyza + gFittyp%lqfit(ka)
      lrt  = gBastyp%kr(la) + gFittyp%nqfit(ka) - 1 - gFittyp%lqfit(ka)
      alfa = gBastyp%alf(la) + gFittyp%alffit(ka)
      nft  = gBasfcn%ndex(lt+2)
      if (nft*nfb>nsinth) call stopit ('overflow SINTh. FA1TTC')

      call mkxyzpointers (lt, 0, 0, ii0, k, icxa)

!     --------------------------------------------------------
!     overlap A (basis*fit) with B (basis): complete cartesian
!     sets
!     --------------------------------------------------------

      if (ka/=llr) then
         if (fopa) then
            call ovrlp (lt, lrt, alfa, lxyzb, lrb, alfb, rab, gBasfcn%ndex, inom, minp, itpow, sint)

!-----------------------------------------------------------------
!**             call fittrform1 (sint, nfb, nft, cxa(icxa), sinth,      
!**                             cxa(icxb))
!-----------------------------------------------------------------

            call mmulrr (sint, nfb, nft, cxa(icxa), nft, nft, 'NO', sinth)
            call mmultr (cxa(icxb), nfb, nfb, sinth, nfb, nft, 'BC', sinth)
         else
            call ovrlp (lxyzb, lrb, alfb, lt, lrt, alfa, rab, gBasfcn%ndex, inom, minp, itpow, sint)

!-----------------------------------------------------------------
!**             call fittrform2 (cxa(icxa), nft, sint, nfb, sinth,      
!**                             cxa(icxb))
!-----------------------------------------------------------------

            call mmultr (cxa(icxa), nft, nft, sint, nft, nfb, 'NO', sinth)
            call mmulrr (sinth, nft, nfb, cxa(icxb), nfb, nfb, 'AC', sinth)
         end if
         llr = ka
      end if

!     --------------------------------------------------------
!     loop over the cartesian fitfunction contributions to fit
!     --------------------------------------------------------

      ctemp(1:nprod) = zero

      nterm_: do nterm = 1, naterm

         id  = numcom(ida)
         idb = kpt + id
         lx  = klist(1,idb)
         ly  = klist(2,idb)
         lz  = klist(3,idb)

         laa = la

!        -------------------------------------------------
!        contribution: select appropriate columns in sinth
!        -------------------------------------------------

         id = 0
         jnfa_: do j = 1, nfa

            lxa = lx + gBastyp%kx(laa)
            lya = ly + gBastyp%ky(laa)
            lza = lz + gBastyp%kz(laa)

            call mkxyzpointers (lxa, lya, lza, ii0, ii, iir)
            ii = ii - ii0
            if (klist(1,ii0+ii)/=lxa) call stopit ('x error')
            if (klist(2,ii0+ii)/=lya) call stopit ('y error')
            if (klist(3,ii0+ii)/=lza) call stopit ('z error')

            if (fopa) then
               ids                = (ii-1)*nfb
               ctemp(id+1:id+nfb) = ctemp(id+1:id+nfb) + cofcom(ida)*sinth(ids+1:ids+nfb)
            else
               iiend             = ii + (nfb-1)*nft
               jend              = j + (nfb-1)*nfa
               ctemp(j:jend:nfa) = ctemp(j:jend:nfa) + cofcom(ida)*sinth(ii:iiend:nft)
            end if

            id = id + nfb
            if (fopa .or. (fopb.and.nft/=1)) laa = laa + 1
         end do jnfa_

         ida = ida + 1
      end do nterm_

!     --------------------------------------------------------------
!     multiply by norms of cartesian basisfunctions and copy to tfit
!     --------------------------------------------------------------

      ndx = iast
      if (fopb) idx = 1

      do ifb = 1, nfb
         fmul = gBastyp%bnorm(lb+ifb-1)
         if (fopa) idx = ifb

         do ifa = 1, nfa
            tfit(ndx) = ctemp(idx)*fmul*gBastyp%bnorm(la+ifa-1)
            ndx       = ndx + nfaa
            idx       = idx + incidx
            
         end do 
      end do 

   end do ia_

!  ---------------------------------------
!  normalize fitfunction and write to file
!  ---------------------------------------

   do iab = 1, nfa*nfb*nfaa, nfaa
      tfit(iab:iab+nfaa-1) = tfit(iab:iab+nfaa-1)*finrm(1:nfaa)
   end do 

   if (fopa) then

      do ifb = 1, nfb
         do ifa = 1,nfa
            k1 = k1 + 1
            indea(k1) = la+ifa - 1
            indeb(k1) = lb+ifb - 1
         enddo
      enddo

   else

      do ifb = 1, nfb
         do ifa = 1,nfa
            k1 = k1 + 1
            do k2 = 1, ncolon
               if (la+ifa-1 .eq. indeb(k2) .and. lb+ifb-1 .eq. indea(k2)) then
                  perm(k1) = k2
               endif
            enddo
         enddo
      enddo

   endif

   k2=1
   do iab = 1, nfa*nfb
      k3 = iab + m0
      if (fopb) k3 = perm(k3)
      do jfaa =  1, nfaa
         amatfit(jfaa+i0, k3) = tfit(k2)
         k2 = k2 + 1
      enddo
   enddo
   m0 = m0 + nfa*nfb

   irec = nfaa*nfa*nfb
   if (irec>gDims%lrl) call stopit ('lrl exceeded. FA1TTC')

!  ----------------
!  end loop over lb
!  ----------------

   lb = lb + nfb
   goto 110

!  ----------------
!  end loop over la
!  ----------------

   210 la = la + nfa
   goto 100

   220 continue
#ifdef D_SCM_MANYTIMERS
   call timere ('fa1ttc')
#endif

end subroutine fa1ttc3c
