! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Main function for replicating Buera and Shin (2013,JPE)
! Author:
!     Xin Tang @ IMF, Spring 2021
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      04/13/2021:          Original Code: No paralleling
!      04/19/2021:            DP Paralleled
!
! Compiling Environment:
!   GNU gfortran on Ubuntu 16.04
!
! Library Used:
!   - MINPACK 
! =============================================================================

program inclusion

    use parameters
    use global
    implicit none

    integer :: indi, indj, inde
    real(dp) :: incr, emean, evar, nummean, numvar
    real(dp) :: wguess, rguess
    real(dp) :: btime, etime
    real(dp), dimension(negrid) :: me, zz, pr
    real(dp), dimension(2) :: x_in, fvec
    integer :: iflag, info, indz, indk, indm
    real(dp), dimension(lwa) :: wa
    external :: fcn_ss_cont

    include 'mpif.h'
    
    ! initialize MPI environment
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    call cpu_time(btime)
    ! Construct grids
    do indi = 1,nintk
        kgrid(indi) = kmin + (indi-1.0_dp)*(kint-kmin)/(nintk-1.0_dp)
    end do

    do indi = nintk+1,nkgrid
        kgrid(indi) = kint + dble(real(indi-nintk))*(kmax-kint) &
         /dble(real(nkgrid-nintk))
    end do

    do indi = 1,nkgrid-1
        do indj = 1,sr+1
            kgridc((indi-1)*sr+indj) = kgrid(indi) + dble(real(indj-1))*&
            (kgrid(indi+1)-kgrid(indi))/dble(real(sr));
        end do
    end do

    do indi = 1,nkgrid-1
        do indj = 1,srsim+1
            kgrids((indi-1)*srsim+indj) = kgrid(indi) + dble(real(indj-1))*&
            (kgrid(indi+1)-kgrid(indi))/dble(real(srsim));
        end do
    end do

    ! entrepreneur productivity
    egrid(negrid)   = exp(-log(1-0.9995)/eta) 
    egrid(negrid-1) = exp(-log(1-0.999)/eta) 
    egrid(negrid-2) = exp(-log(1-0.998)/eta)
    egrid(1)        = exp(-log(1-0.3670)/eta)
    incr = (egrid(negrid-2)-egrid(1))/dble(real(negrid-2-1))
    do inde = 2,negrid-2
        egrid(inde) = egrid(inde-1)+incr 
    end do
    me = 1-egrid**(-eta) 
    do inde = 2,negrid
        edist(inde) = (me(inde) - me(inde-1))/me(negrid)
    end do
    edist(1) = me(1)/me(negrid)

    emean = eta/(eta-1)
    evar = (1/(eta-1))**2*(eta/(eta-2)) 
    nummean = sum(edist*egrid)
    numvar = sum(((egrid-nummean)**2)*edist)

    egrid = egrid*0.2_dp

    omega(1:negrid) = edist*(1-exp(-q*egrid)) ;
    omega(negrid+1:nee) = edist*(exp(-q*egrid)) ;
    tauz(:,1) = (1-tplus)*egrid ;
    tauz(:,2) = (1-tminus)*egrid ;

    if (myrank .eq. root) then 
        write (*,*) 'Pareto distribution: compare mean and variance'
        write (*,*) 'Theoretical: ', emean, evar
        write (*,*) 'Numerical', nummean, numvar
    end if

    ! Distribute workload
    indm = 0
    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            indm = indm + 1
            iifun(indk,(indz-1)*negrid+inde) = indm
            izfun(indm) = indz
            iefun(indm) = inde
            ikfun(indm) = indk
        end do
    end do
    end do

    ni_indi = int(real(ni_total-1)/real(nproc))+1
    allocate(mvali(ni_indi),mpolci(ni_indi),mpolai(ni_indi),mpolaindi(ni_indi))
    allocate(mvala(ni_total),mpolca(ni_total),mpolaa(ni_total),mpolainda(ni_total))
    ibegin = myrank*ni_indi + 1
    iend = min((myrank+1)*ni_indi,ni_total)

    wguess = 0.0803852580114_dp
    rguess = -0.0501250_dp
    x_in(1) = wguess
    x_in(2) = rguess

    ! evaluate excess demand at price xin
    ! call fcn_ss(x_in,fvec)
    ! call fcn_ss_cont(nmarket,x_in,fvec,iflag)
    call hybrd1(fcn_ss_cont,nmarket,x_in,fvec,tolmin,info,wa,lwa)

    if (myrank .eq. root) then
        !write (*,*) 'wage == ', wguess, 'interest == ', rguess
        write (*,*) 'wage == ', x_in(1), 'interest == ', x_in(2)
        write (*,*) 'excess demand'
        write (*,*) 'labor == ', fvec(1), 'capital == ', fvec(2)
    end if

    if (myrank .eq. root) then
        call save_results
    end if

    call cpu_time(etime)

    if (myrank .eq. root) then
        write (*,*) 'time running == ', etime-btime, ' secs'
    end if

    deallocate(mvali,mpolci,mpolai,mpolaindi)
    deallocate(mvala,mpolca,mpolaa,mpolainda)

    call mpi_finalize(ierr)

end program inclusion