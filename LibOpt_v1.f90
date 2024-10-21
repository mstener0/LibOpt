program OpenLib_Ak

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This program calculates the matrix Ak for the algorithm PolTDDFT:                      !
    !                                                                                          !
    !   https://www.youtube.com/watch?v=Y0AvVS1bQKM   video itroduction by Profesor M. Stener  !
    !   https://pubs.aip.org/aip/jcp/article/143/2/024106/825263              related article  !
    !                                                                                          !
    !   The expression for the generic term of Ak is:                                          ! 
    !                                                                                          !
    !   Ak(mu,ia) = < mu | ia >                                                                !
    !             = Sum_over_nbas_index_sigma(Sum_over_nbas_index_tau(< mu | sigma*tau > *     !
    !               * C_sigma_i * C_tau_a))                                                    !
    !                                                                                          !
    !   where:                                                                                 !
    !          mu        is an element of the density fitting functions basis set              !
    !          i         is an occupied molecular orbital                                      !
    !          a         is a virtual molecular orbital                                        !
    !          nbas      is the number of atomic orbital basis functions                       !
    !          sigma     is an element of the atomic orbital basis functions set used to       !
    !                    describe occupied molecular orbital i                                 !
    !          tau       is an element of the atomic orbital basis functions set used to       !
    !                    describe occupied molecular orbital a                                 !
    !          C_sigma_i it is the weight of atomic orbital sigma in the description of the    !
    !                    occupied molecular orbital i                                          !
    !          C_tau_a   it is the weight of atomic orbital tau in the description of the      !
    !                    virtual molecular orbital a                                           !
    !                                                                                          !
    !   The calculation of such matrix is computationally demanding and to speed up the proces !
    !   the Ak matrix will be calculated in parallel, using GPUs, by dividing it into multiple !
    !   submatrices, one for each fitting function mu.                                         !
    !   The main objects needed for the calculation are the integrals < mu | sigma*tau > and   !
    !   the eigenvalues C_sigma_i, C_tau_a. Such objects are given as input and used to build  !
    !   matrices mu_pairfit(/final_dim_mu_pairfit, /final_dim_mu_pairfit),                     !
    !   mu_eigocc(/final_dim_mu_pairfit, /mx_occ) and eigvirt(/final_dim_mu_pairfit, /mx_virt) !
    !   Where:                                                                                 !
    !          mu_pairfit   is a matrix of reals with the integrals beetwen the specific mu    !
    !                       fit functions and the atomic orbital basis functions sigma*tau.    !
    !                       The matrix is square having along the rows sigma and columns tau.  !
    !                       Only bas functions able to create at least one couple interacting  !
    !                       with the specific mu fitting functions are considered              !
    !                                                                                          !
    !          mu_eigocc    is a matrix of reals with the coefficients of the expansion of the !
    !                       occupied molecular orbitals (same index of the columns) over the   !
    !                       atomic orbitals basis functions (but only basis functions able to  !
    !                       create at least one couple interacting with the specific mu fit    !
    !                       function are considered).                                          !
    !                                                                                          !
    !          mu_eigvirt   same as above but for virtual molecular orbitals                   !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer :: nmu, mx_ncpl, nbas, mx_occ, mx_virt, io_status
      ! nmu is the number of fitting functions "mu", it is read from input file
      ! mx_ncpl is the maximum number of bas function couples among all fit functions, it is
      !          read from input and used to allocate matrices with enougth space
      ! nbas is the number of atomic orbitals basis functions from witch we pick "sigma" and 
      !       "tau", it is read from input file
      ! mx_occ is the maximum occupied molecular orbital, it is read from input file
      ! mx_virt is the maximum virtual molecular orbital, it is read from input file

      integer :: mu_mx_ncpl, mu_mx_basa, mu_mx_basb, mu_nbasa, mu_nbasb
      ! mu_mx_ncpl is the number of bas function couples for a specific mu fitting function,
      !            it is calculated in "mx_for_mu" subroutine
      ! mu_mx_basa is the maximum abs. index of bas function "sigma" for the bas couples able
      !            to interact with specific fit function mu, it is calculated in "mx_for_mu"
      !            subroutine
      ! mu_mx_basb same as above but for bas function "tau"
      ! mu_nbasa it is the number of different "sigma" functions considered in the couples of
      !          bas functions able to interact with specific fit function mu, it is
      !          calculated in "mx_for_mu" subroutine
      ! mu_nbasb same as above but for bas function "tau"
      
      integer :: mu_mx_bas, i, j, final_dim_mu_pairfit, mu
      ! mu_mx_bas is the maximum value beetwen mu_mx_basa and mu_mx_bas
      ! final_dim_mu_pairfit is the dimennsion of squared matrix mu_pairfit_tmp, it is 
      !                      calculated in "final_dim_mu_pairfit" subroutine
      ! mu is an index to run over the fitting functions
      ! i & j are multiporpouse indeces to run over different loops

      integer, allocatable :: basa(:,:), basb(:,:)
      ! basa is a matrix that contains on its columns the abs. indeces of atomic orbital bas
      !      functions "sigma" involved in the bas couples able to interact with fit function
      !      mu (having index the same row index of such matrix). It is created in subroutine
      !      "read_matrices_data" using data from input file
      ! basb same as above but for bas function "tau"

      real*8, allocatable :: pairfit(:,:), eigocc(:,:), eigvirt(:,:)
      ! pairfit is a matrix that contains the value of integrals < mu | sigma*tau > with 
      !         index of mu equal to the index of the row of such matrix and along the
      !         columns are all the possible bas function couples able to interact with mu
      !         fit function. It is created in subroutine "read_matrices_data" using data
      !         from input file
      ! eigocc is a matrix with the coefficients of the expansion of the occupied molecular
      !        orbitals (the columns) over the set of atomic orbitals basis functions (the
      !        rows). It is created in subroutine "read_matrices_data" using data from input
      !        file
      ! eigvirt same as above but for virtual molecular orbitals

      real*8, allocatable :: mu_pairfit_tmp(:,:), mu_pairfit(:,:), mu_eigocc(:,:),&
                           & mu_eigvirt(:,:), mu_pfeigocc(:,:), mu_eigvirt_trp(:,:),&
                           & mu_Ak(:,:), mu_Ak_trp(:,:) 
                           
      ! mu_pairfit_tmp is a squared matrix specific for each fit function mu thet reports the
      !            integrals < mu | sigma*tau > accordingly with (abs. index sigma) = index
      !            of the row and (abs. index tau) = index of the column. This way many rows
      !            and columns are empty. It is created in "build_pairfit_for_fit" subroutine
      ! mu_pairfit it is mu_pairfit_tmp cleaned from empty rows and columns. It is created
      !            in "final_pairfit_for_fit" subroutine
      ! mu_eigocc it is a matrix like eigocc but specific for the fit function mu from witch
      !           the rows related to bas functions that do not appears in any bas function
      !           couple able to interact with specific mu fit function are eliminated. It is
      !           created in "mu_eigen" subroutine
      ! mu_eigvirt same as above but for eigvirt
      ! mu_pfeigocc is the matrix product of mu_pairfit*mu_eigocc
      ! mu_eigvirt_trp is the transposed matrix of mu_eigvirt
      ! mu_Ak is the submatrix Ak specific for this mu fit function, on the rows are the
      !          virtual molecular orbitals and on the columns the occupied ones

      real*8, allocatable :: bas_to_trim(:)
      ! bas_to_trim is a vector used to tell if a bas function has to be considered or not
      !             for a specific fit function mu. It is created in "build_pairfit_for_fit"
      !             subroutine

      character(len=256) :: filename
      ! filename is the path/name of the input file


      if (command_argument_count() < 1) then
         print *, "Usage: ./OpenLib_Ak <path_to_input_file>"
         stop
      endif

      call get_command_argument(1, filename)

      call read_matrices_dimensions(filename, nmu, mx_ncpl, nbas, mx_occ, mx_virt)
      
      allocate(basa(nmu, mx_ncpl))
      allocate(basb(nmu, mx_ncpl))
      allocate(pairfit(nmu, mx_ncpl))
      allocate(eigocc(nbas, mx_occ))
      allocate(eigvirt(nbas, mx_virt))

      call read_matrices_data(filename, nmu, mx_ncpl, nbas, mx_occ, mx_virt, basa,&
                             & basb, pairfit, eigocc, eigvirt)
      OPEN(UNIT=10, FILE="Ak_matrices.txt", STATUS='REPLACE', ACTION='WRITE', IOSTAT=io_status)

      if (io_status /= 0) THEN
         PRINT *, 'Error opening file Ak_matrices.txt'
         STOP
      endif

      do mu = 1, nmu
         call mx_for_mu(mx_ncpl, basa(mu,:), basb(mu,:), mu_mx_ncpl, mu_mx_basa,&
                       & mu_mx_basb, mu_nbasa, mu_nbasb)
         mu_mx_bas = max(mu_mx_basa,mu_mx_basb)
         
         allocate(mu_pairfit_tmp(mu_mx_bas, mu_mx_bas))

         mu_pairfit_tmp = 0.d0

         allocate(bas_to_trim(mu_mx_bas))

         bas_to_trim = 0.d0
         call build_pairfit_for_fit(nbas, mu_mx_ncpl, basa(mu,1:mu_mx_ncpl),&
                                   & basb(mu,1:mu_mx_ncpl),&
                                   & pairfit(mu,1:mu_mx_ncpl), mu_mx_bas, mu_pairfit_tmp,&
                                   & final_dim_mu_pairfit, bas_to_trim)
         allocate(mu_pairfit(final_dim_mu_pairfit, final_dim_mu_pairfit))

         mu_pairfit = 0.d0

         call final_pairfit_for_fit(nbas, mu_pairfit_tmp, mu_mx_bas, mu_pairfit,&
                                   & final_dim_mu_pairfit, bas_to_trim)
        
         allocate(mu_eigocc(final_dim_mu_pairfit,mx_occ))
         allocate(mu_eigvirt(final_dim_mu_pairfit,mx_virt))
         allocate(mu_eigvirt_trp(mx_virt, final_dim_mu_pairfit))
         allocate(mu_Ak(mx_virt, mx_occ))

         mu_eigocc = 0.d0
         mu_eigvirt = 0.d0
         mu_eigvirt_trp = 0.d0
         mu_Ak = 0.d0

         call mu_eigen(eigocc, eigvirt, mx_occ, mx_virt, nbas, bas_to_trim,&
                      & final_dim_mu_pairfit, mu_mx_bas, mu_eigocc, mu_eigvirt)

         allocate(mu_pfeigocc(final_dim_mu_pairfit, mx_occ))

         mu_pfeigocc = 0.d0

         call mu_pairfit_x_eigocc(mu_pairfit, mu_eigocc, final_dim_mu_pairfit, mx_occ,&
                                 & mu_pfeigocc)

         do i = 1, final_dim_mu_pairfit
            do j = 1, mx_virt
               mu_eigvirt_trp(j, i) = mu_eigvirt (i, j)
            enddo
         enddo

         call mu_eigvirt_trp_x_mu_pfeigocc(mu_pfeigocc, mu_eigvirt_trp,&
                                          & final_dim_mu_pairfit, mx_virt, mx_occ, mu_Ak)

         allocate(mu_Ak_trp(mx_occ, mx_virt))

         mu_Ak_trp = 0.d0

         do i = 1, mx_virt
            do j = 1, mx_occ
               mu_Ak_trp(j,i) = mu_Ak(i,j)
            enddo
         enddo

         write(10,*) "                                 "
         write(10,*) "Ak -", mu
         write(10,*) "                                 "
         do i = 1, mx_occ
            do j = 1, mx_virt
               write(10,*) "Ak-",mu,"(",i,",",j,") = ", mu_Ak_trp(i,j)
            enddo
         enddo

         deallocate(mu_pairfit_tmp)
         deallocate(mu_pairfit)
         deallocate(bas_to_trim)
         deallocate(mu_eigocc)
         deallocate(mu_eigvirt)
         deallocate(mu_pfeigocc)
         deallocate(mu_eigvirt_trp)
         deallocate(mu_Ak)
         deallocate(mu_Ak_trp)
      enddo

      CLOSE(UNIT=10)

      deallocate(basa)
      deallocate(basb)
      deallocate(pairfit)
      deallocate(eigocc)
      deallocate(eigvirt)

end program OpenLib_Ak

subroutine read_matrices_dimensions(filename, nmu, mx_ncpl, nbas, mx_occ, mx_virt)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine reads from filename (given as an input when running the code):         !
    !                                                                                          !
    !   nmu: the number of fitting functions mu                                                !
    !   mx_ncpl: the maximum number of basis function couples i.e. sigma*tau among all mu      !
    !   nbas: the number of basis functions                                                    !
    !   mx_occ: the index of the maximum occupied molecular orbital                            !
    !   mx_virt: the index of the maximum virtual molecular orbital                            !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    
    integer, intent(out) :: nmu, mx_ncpl, nbas, mx_occ, mx_virt
    integer :: unit_num
    character(len=*), intent(in) :: filename

    unit_num = 20
    open(unit=unit_num, file=filename, status='old', action='read')

    read(unit_num, *) nmu, mx_ncpl, nbas, mx_occ, mx_virt

    close(unit_num)

end subroutine read_matrices_dimensions

subroutine read_matrices_data(filename, nmu, mx_ncpl, nbas, mx_occ, mx_virt, basa, basb,&
                             & pairfit, eigocc, eigvirt)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine builds from filename the matrices:                                     !
    !                                                                                          !
    !   basa(/nmu, /mx_ncpl): it is a matrix of integers. It reports, for every fitting        !
    !                         function mu, the abs. index of the bas function sigma            !
    !                                                                                          !
    !   basb(/nmu, /mx_ncpl): same as above but for the bas function tau                       !
    !                                                                                          !
    !   pairfit(/nmu, /mx_ncpl): for every fitting function mu it reports, for each bas        !
    !                            function couple sigma*tau, the real value of the integral     !
    !                            < mu | sigma*tau >                                            !
    !                                                                                          !
    !   All the above matrices has some 0s at the end of each row becase that dimension is     !
    !   allocated at the number of bas couples possesed by the fit function mu having the most !
    !   of them.                                                                               ! 
    !                                                                                          !
    !   eigocc(/nbas, /mx_occ): it is a matrix of reals. It reports for every bas function its !
    !                           contribution to the description of each occupied molecular     !
    !                           orbital.                                                       !
    !                                                                                          !
    !   eigvirt(/nbas, /mx_virt): same as before but for virtual molecular orbitals            !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    
    integer, intent(in) :: nmu, mx_ncpl, nbas, mx_occ, mx_virt
    integer, intent(out) :: basa(nmu, mx_ncpl), basb(nmu, mx_ncpl)
    real*8, intent(out) :: pairfit(nmu, mx_ncpl), eigocc(nbas, mx_occ), eigvirt(nbas, mx_virt)
    integer :: i, j
    integer :: unit_num
    character(len=256) :: line 
    character(len=*), intent(in) :: filename

    unit_num = 20
    open(unit=unit_num, file=filename, status='old', action='read')

    read(unit_num, *)

    read(unit_num, '(A)') line

    do i = 1, nmu
        read(unit_num, '(10I6)') (basa(i, j), j = 1, mx_ncpl)
    end do

    read(unit_num, '(A)') line

    do i = 1, nmu
        read(unit_num, '(10I6)') (basb(i, j), j = 1, mx_ncpl)
    end do

    read(unit_num, '(A)') line

    do i = 1, nmu
        read(unit_num, '(5ES23.15E3)') (pairfit(i, j), j = 1, mx_ncpl)
    end do

    read(unit_num, '(A)') line

    do i = 1, nbas
        read(unit_num, '(5ES23.15E3)') (eigocc(i, j), j = 1, mx_occ)
    end do

    read(unit_num, '(A)') line

    do i = 1, nbas
        read(unit_num, '(5ES23.15E3)') (eigvirt(i, j), j = 1, mx_virt)
    end do

    close(unit_num)

end subroutine read_matrices_data

subroutine mx_for_mu(mx_ncpl, mu_basa, mu_basb, mu_mx_ncpl, mx_basa, mx_basb, nbasa, nbasb)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine calculates: mu_mx_ncpl (the number of basis function couples for fit   !
    !   function mu), mu_nbasa and mu_nbasb (the number of different basis functions sigma and !
    !   tau involved in a couple able to interact with fit function mu), mx_basa and mx_basb   !
    !   are the maximum abs. index of such functions sigma and tau respectively.               !
    !   They are calculated running over mx_ncpl and accumulating on mu_mx_ncpl                !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer, intent(in)  :: mx_ncpl
    integer, intent(in)  :: mu_basa(mx_ncpl), mu_basb(mx_ncpl)
    integer, intent(out) :: nbasa, nbasb
    integer, intent(out) :: mu_mx_ncpl, mx_basa, mx_basb
    integer              :: supporting_basa(mx_ncpl), supporting_basb(mx_ncpl)
    integer              :: i

    mu_mx_ncpl = 0
    supporting_basa = 0
    supporting_basb = 0
    mx_basa = 0
    mx_basb = 0

    do i = 1, mx_ncpl
       if (mu_basa(i).ne.0) then
          mu_mx_ncpl = mu_mx_ncpl + 1
          supporting_basa(mu_basa(i)) = 1
          supporting_basb(mu_basb(i)) = 1
          if (mu_basa(i).gt.mx_basa) mx_basa = mu_basa(i)
          if (mu_basb(i).gt.mx_basb) mx_basb = mu_basb(i)
       else
          exit
       endif
    enddo

    nbasa = SUM(supporting_basa)
    nbasb = SUM(supporting_basb)

end subroutine mx_for_mu

subroutine build_pairfit_for_fit(nbas, mu_mx_ncpl, basa, basb, pairfit, mu_mx_bas,&
                                & mu_pairfit_tmp, final_dim_mu_pairfit, bas_to_trim)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine builds mu_pairfit_tmp, bas_to_trim and calculates final_dim_mu_pairfit !
    !                                                                                          !
    !   mu_pairfit_tmp: it is a square matrix with dimension mu_mx_bas. It is created rugging  ! 
    !                   over the couples of basis functions that interact with fit function mu !
    !                   (there are mu_mx_ncpl of them) and using vectors basa and basb to know !
    !                   abs. index of the basis functions involved in such a couple.           !
    !                                                                                          !
    !   bas_to_trim: it is a vector of reals with dimension mu_mx_bas. It is created by runnig !
    !                over rows and columns of mu_pairfit_tmp to check what couples of basis    !
    !                functions do not interact with fit function mu. This way bas_to_trim has  !
    !                value 0.d0 if that row and column of mu_pairfit_tmp can be eliminated, it !
    !                has value 1.d0 otherwise.                                                 !
    !                                                                                          !
    !   final_dim_mu_pairfit: is the dimension of final_dim_mu_pairfit (mu_pairfit_tmp cleaned !
    !                         of empty rows and columns). It is created adding toghether the   !
    !                         elements of bas_to_trim.                                         !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer, intent(in) :: mu_mx_ncpl, nbas, mu_mx_bas
    integer, intent(in) :: basa(mu_mx_ncpl), basb(mu_mx_ncpl)
    integer, intent(out):: final_dim_mu_pairfit
    real*8, intent(in)  :: pairfit(mu_mx_ncpl)
    real*8, intent(out) :: mu_pairfit_tmp(mu_mx_bas, mu_mx_bas)   
    real*8, intent(out) :: bas_to_trim(mu_mx_bas)
    integer :: i_cpl, mx_basa, mx_basb, i

    mu_pairfit_tmp = 0.d0
    bas_to_trim= 0.d0

    do i_cpl = 1, mu_mx_ncpl
       mu_pairfit_tmp(basa(i_cpl),basb(i_cpl)) = pairfit(i_cpl)
       mu_pairfit_tmp(basb(i_cpl),basa(i_cpl)) = pairfit(i_cpl)
    enddo

    do i=1,mu_mx_bas
       bas_to_trim(i)= SUM(mu_pairfit_tmp(:,i)) + SUM(mu_pairfit_tmp(i,:))
       if (bas_to_trim(i).ne.0.d0) bas_to_trim(i)=1
    enddo
    
    final_dim_mu_pairfit = INT(SUM(bas_to_trim))

end subroutine build_pairfit_for_fit
       
subroutine final_pairfit_for_fit(nbas, mu_pairfit_tmp, mu_mx_bas, mu_pairfit,&
                                & final_dim_mu_pairfit, bas_to_trim)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine builds mu_pairfit by removing from mu_pairfit_tmp the rows and         !
    !   columns related to basis functions that do not interact with fitting function mu       !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer, intent(in) :: mu_mx_bas, final_dim_mu_pairfit, nbas
    real*8, intent(in)  :: mu_pairfit_tmp(mu_mx_bas, mu_mx_bas), bas_to_trim(mu_mx_bas)
    real*8, intent(out) :: mu_pairfit(final_dim_mu_pairfit,final_dim_mu_pairfit)
    real*8   :: tmp_a_mu_pairfit_tmp(final_dim_mu_pairfit, mu_mx_bas)
    integer  :: i, l

    l = 0

    do i = 1, mu_mx_bas
       if (bas_to_trim(i).ne.0.d0) then
          l = l + 1
          tmp_a_mu_pairfit_tmp(l,:) = mu_pairfit_tmp(i,:)
       endif
    enddo

    l = 0

    do i = 1, mu_mx_bas
       if (bas_to_trim(i).ne.0.d0) then
          l = l + 1
          mu_pairfit(:,l) = tmp_a_mu_pairfit_tmp(:,i)
       endif
    enddo

end subroutine final_pairfit_for_fit

subroutine mu_eigen(eigocc, eigvirt, mx_occ, mx_virt, nbas, bas_to_trim, final_dim_mu_pairfit,&
                   & mu_mx_bas, mu_eigocc, mu_eigvirt)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine builds mu_eigocc and mu_eigvirt by removing from eigocc and eigvirt    !
    !   the rows related to basis funtions that do not interact with fitting function mu       !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer, intent(in) :: mx_occ, mx_virt, nbas, final_dim_mu_pairfit, mu_mx_bas
    real*8, intent(in)  :: eigocc(nbas, mx_occ), eigvirt(nbas, mx_virt), bas_to_trim(mu_mx_bas)
    real*8, intent(out) :: mu_eigocc(final_dim_mu_pairfit, mx_occ),&
                           & mu_eigvirt(final_dim_mu_pairfit, mx_virt)
    integer  :: l, i

    l = 0

    mu_eigocc = 0.d0
    mu_eigvirt = 0.d0

    do i = 1, mu_mx_bas
       if (bas_to_trim(i).ne.0.d0) then
          l = l + 1
          mu_eigocc(l,:) = eigocc(i,:)
          mu_eigvirt(l,:) = eigvirt(i,:)
       endif
    enddo

end subroutine mu_eigen

subroutine mu_pairfit_x_eigocc(mu_pairfit, mu_eigocc, final_dim_mu_pairfit, mx_occ,&
                              & mu_pfeigocc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine execute the matrix multiplication mu_pairfit*mu_eigocc and returns the !
    !   resulting matrix as mu_pfeigocc                                                        !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer, intent(in) :: final_dim_mu_pairfit, mx_occ
    real*8, intent(in)  :: mu_pairfit(final_dim_mu_pairfit, final_dim_mu_pairfit),&
                           mu_eigocc(final_dim_mu_pairfit, mx_occ)
    real*8, intent(out) :: mu_pfeigocc(final_dim_mu_pairfit, mx_occ)

    integer  :: i, j, k

    do i = 1, final_dim_mu_pairfit
       do j = 1, mx_occ
          do k = 1, final_dim_mu_pairfit
             mu_pfeigocc(i,j) = mu_pfeigocc(i,j) + mu_pairfit(i,k) * mu_eigocc(k,j)
          end do
       end do
    end do

end subroutine mu_pairfit_x_eigocc

subroutine mu_eigvirt_trp_x_mu_pfeigocc(mu_pfeigocc, mu_eigvirt_trp, final_dim_mu_pairfit,&
                                       & mx_virt, mx_occ, mu_Ak)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          !
    !   This subroutine execute the matrix multiplication mu_eigvirt_trp*mu_pfeigocc, the      !
    !   resulting matrix is mu_Ak, the specific submatrix of Ak related to fit function mu     !                                                      !
    !                                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer, intent(in) :: final_dim_mu_pairfit, mx_occ, mx_virt
    real*8, intent(in)  :: mu_pfeigocc(final_dim_mu_pairfit, mx_occ),&
                           & mu_eigvirt_trp(mx_virt, final_dim_mu_pairfit)
    real*8, intent(out) :: mu_Ak(mx_virt, mx_occ)

    integer  :: i, j, k

    do i = 1, mx_virt
       do j = 1, mx_occ
          do k = 1, final_dim_mu_pairfit
                mu_Ak(i,j) = mu_Ak(i,j) + mu_eigvirt_trp(i,k) * mu_pfeigocc(k,j)
          end do
       end do
    end do

end subroutine mu_eigvirt_trp_x_mu_pfeigocc
