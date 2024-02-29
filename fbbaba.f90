#include <config.h>
subroutine fbbaba (nfaa, lalo, lahi, kalo, kahi, iskf, numcom, cofcom, rab,       &
                       inom, minp, itpow, fopa, ncxa, cxa,    &
                       nalo, irec, tfitfbbaba)

   use KF
   use ConstModule
   use FactorialsModule
   use MasterModule
   use DimensionsModule
!   use FitIntBuffer


   use Vartypes
   implicit real(KREAL) (a-h, o-z)
   implicit integer(KINT) (i-n)


   integer(KINT), intent(in)  :: nfaa, lalo, lahi, kalo, kahi, nalo, irec
   real(KREAL), intent(in)    :: rab

   integer(KINT), intent(in)  :: inom(496), minp(31), itpow(31), iskf(4,gDims%niskf),              &
          numcom(gDims%na1cof)
   logical, intent(in)        :: fopa
   real(KREAL), intent(in)    :: cxa(ncxa), cofcom(gDims%na1cof)
   real(KREAL), intent(inout) :: tfitfbbaba(irec)

!  ======================================================================
!  purpose: calculate the two center fit integrals using
!           A1 fit functions
!
!  input: nfaa - number of fitfunctions on center A
!         lalo - lower index for basis functions on center A
!         lahi - upper index for basis functions on center A
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
!         nalo - lower index in sets of basis functions on a
!
!  remark: the calculated fit integrals are written to the current
!          variable of unit iuf as packed reals, thus possibly extending
!          that variable
!
!  called from: fitint
!  ======================================================================

   real(KREAL), parameter   :: zero = 0.0_KREAL
   integer(KINT), parameter :: nctemp = 150, nsinth = 420, nklist = 120, lklist = 7

   integer(KINT) :: klist(3,nklist)
   logical       :: fopb, firstk
   real(KREAL)   :: ctemp(nctemp), sinth(nsinth), sint(nsinth) !prova 

   real(KREAL), allocatable :: tfit(:) !prova
   integer(KINT) :: max_nfa, max_nfb !prova 

   save firstk,klist

   data firstk /.true./

#ifdef D_SCM_MANYTIMERS
   call timers ('fbbaba')
#endif


   if (firstk) then
      firstk = .false.
      call mkxyzlist (lklist, nklist, k, klist)
      if (k/=nklist) call stopit ('klist error. FA1TTC')
   end if

   fopb   = .not.fopa
   incidx = 1

   max_nfa = 0_KINT
   max_nfb = 0_KINT

!  ------------------------------------------
!  loop over the sets of basis functions on A
!  ------------------------------------------

   iseta = nalo - 1
   la    = lalo
300 if (la>lahi) goto 420
   iseta = iseta + 1

   lxyza = gBastyp%kx(la) + gBastyp%ky(la) + gBastyp%kz(la)
   nfa   = gBasfcn%ndex(lxyza+2)
   max_nfa = max(max_nfa, nfa)

!  ------------------------------------------
!  loop over the sets of basis functions on A
!  ------------------------------------------

   isetb = nalo - 1
   lb    = lalo
310 if (lb>la) goto 410
   isetb = isetb + 1

   lxyzb = gBastyp%kx(lb) + gBastyp%ky(lb) + gBastyp%kz(lb)
   nfb   = gBasfcn%ndex(lxyzb+2)
   max_nfb = max(max_nfb, nfb)
   lrb   = gBastyp%kr(lb)
   alfb  = gBastyp%alf(lb)

!  ----------------
!  end loop over lb
!  ----------------

   lb = lb + nfb

   goto 310

!  ----------------
!  end loop over la
!  ----------------

410 la = la + nfa

   goto 300

420 continue

   allocate(tfit(nfaa*max_nfa*max_nfb))

!  ------------------------------------------
!  loop over the sets of basis functions on A
!  ------------------------------------------

   iseta = nalo - 1
   la    = lalo
100 if (la>lahi) goto 220
   iseta = iseta + 1

   lxyza = gBastyp%kx(la) + gBastyp%ky(la) + gBastyp%kz(la)
   nfa   = gBasfcn%ndex(lxyza+2)

!  ------------------------------------------
!  loop over the sets of basis functions on A
!  ------------------------------------------

   isetb = nalo - 1
   lb    = lalo
110 if (lb>la) goto 210
   isetb = isetb + 1

   lxyzb = gBastyp%kx(lb) + gBastyp%ky(lb) + gBastyp%kz(lb)
   nfb   = gBasfcn%ndex(lxyzb+2)
   lrb   = gBastyp%kr(lb)
   alfb  = gBastyp%alf(lb)

!  ------------------------------------
!  make sure the buffer is large enough
!  ------------------------------------

   if (fopa) incidx = nfb

   llr   = -1
   nprod = nfa*nfb
   if (nprod>nctemp) call stopit ('overflow CTEMP. FA1TTC')

!  --------------------------------------------------------
!  combine basis*basis function characteristics
!  --------------------------------------------------------

   lt   = lxyza + lxyzb  !  gFittyp%lqfit(ka)
   lrt  = gBastyp%kr(la) + gBastyp%kr(lb) ! gFittyp%nqfit(ka) - 1 - gFittyp%lqfit(ka)
   alfa = gBastyp%alf(la) + gBastyp%alf(lb)  ! gFittyp%alffit(ka)
   nft  = gBasfcn%ndex(lt+2)

   call mkxyzpointers (lt, 0, 0, ii0, k, icxa)

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
      nff = gBasfcn%ndex(gFittyp%lqfit(ka)+2)
      if (nft*nff>nsinth) call stopit ('overflow SINTh. FA1TTC')
      call mkxyzpointers (gFittyp%lqfit(ka), 0, 0, k0, k, icxb)

!     --------------------------------------------------------
!     overlap A (basis*basis) with B (fit): complete cartesian
!     sets
!     --------------------------------------------------------

      if (ka/=llr) then
         call ovrlp (lt, lrt, alfa, gFittyp%lqfit(ka), gFittyp%nqfit(ka) - 1 - gFittyp%lqfit(ka),  gFittyp%alffit(ka), &
                    & rab, gBasfcn%ndex, inom, minp, itpow, sint)

         call mmulrr (sint, nff, nft, cxa(icxa), nft, nft, 'NO', sinth)
         call mmultr (cxa(icxb), nff, nff, sinth, nff, nft, 'BC', sinth)

         llr = ka
      end if

!     --------------------------------------------------------
!     loop over the cartesian fitfunction contributions to fit
!     --------------------------------------------------------

      ctemp(1:nprod) = zero

      do i = 1, nft
         do j = 1, naterm
            idx = (i-1)*nff+numcom(ida-1+j)

            ctemp(i) = ctemp(i) + cofcom(ida-1+j)*sinth(idx)

         enddo
      enddo

      do i=1,nfa
         laa = la -1 + i
         result1 = gBastyp%bnorm(laa)
         do j=1,nfb
            lbb = lb - 1 + j
            if (lbb>laa) exit

            lxa = gBastyp%kx(laa) + gBastyp%kx(lbb)
            lya = gBastyp%ky(laa) + gBastyp%ky(lbb)
            lza = gBastyp%kz(laa) + gBastyp%kz(lbb)

            call mkxyzpointers (lxa, lya, lza, ii0, ii, iir)
            ii = ii - ii0

            if (klist(1,ii0+ii)/=lxa) call stopit ('x error')
            if (klist(2,ii0+ii)/=lya) call stopit ('y error')
            if (klist(3,ii0+ii)/=lza) call stopit ('z error')

            result = ctemp(ii) * gBastyp%bnorm(lbb) * result1

            ndx = (i-1)*(nfb*nfaa)+(j-1)*nfaa+iast
            tfit(ndx)=result
            ndxb = ndx

         enddo
      enddo

   end do ia_

   do ifb = lb, lbb
      do ifa = la, laa 
     
         if (ifb > ifa) cycle

         infa = ifa - lalo + 1 
         infb = ifb - lalo + 1

         indice = ((ifa-la)*(lbb-lb+1)+(ifb-lb))*nfaa +1

         in = 1 + nfaa*(infa*(infa-1)/2 + infb - 1)
         ifin = nfaa*(infa*(infa-1)/2 + infb)

         tfitfbbaba( in : ifin) = tfit(indice:indice-1+nfaa)     

      enddo
   enddo 

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

   deallocate(tfit)

#ifdef D_SCM_MANYTIMERS
   call timers ('fbbaba')
#endif
end subroutine fbbaba

