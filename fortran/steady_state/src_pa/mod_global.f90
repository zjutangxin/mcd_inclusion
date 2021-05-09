module global

    use parameters
    implicit none

    real(dp), dimension(nkgrid) :: kgrid
    real(dp), dimension(nkgridc) :: kgridc
    real(dp), dimension(nsim) :: kgrids
    real(dp), dimension(negrid) :: egrid, edist
    real(dp), dimension(nee) :: omega
    real(dp), dimension(negrid,2) :: tauz

    real(dp), dimension(nkgrid,nee) :: polc, pola, polk, poll, poly, polpi, vn
    real(dp), dimension(nsim,nee) :: dpolc, dpola, dpolk, dpoll, dpoly, dpolpi
    integer, dimension(nsim,nee) :: dpolo
    integer, dimension(nkgrid,nee) :: polaind, polo
    real(dp), dimension(nsim,nee) :: ddistss
    real(dp) :: ademand, asupply, ldemand, lsupply, wss, rss
    real(dp), dimension(nee) :: wgt, maxprof

    ! aggregate statistics
    real(dp) :: govtax, govsub, govex, rcex
    real(dp) :: aggy, aggc, aggk, aggef, aggl, nfirm, emeanf
    real(dp) :: lmean, yperfirm, lprod, pimean, incmean
    real(dp) :: share_ld, share_id, fexit

    ! MPI variables
    integer :: ierr, myrank, nproc, ni_indi, ibegin, iend
    real(dp), dimension(:), allocatable :: mvali, mvala
    real(dp), dimension(:), allocatable :: mpolci, mpolai, mpolca, mpolaa
    integer, dimension(:), allocatable :: mpolaindi, mpolainda
    integer, dimension(ni_total) :: izfun, iefun, ikfun
    integer, dimension(nkgrid,nee) :: iifun

end module global