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
    real(dp) :: ademand, asupply, ldemand, lsupply

    ! MPI variables
    integer :: ierr, myrank, nproc, ni_indi, ibegin, iend
    real(dp), dimension(:), allocatable :: mvali, mvala
    real(dp), dimension(:), allocatable :: mpolci, mpolai, mpolca, mpolaa
    integer, dimension(:), allocatable :: mpolaindi, mpolainda
    integer, dimension(ni_total) :: izfun, iefun, ikfun
    integer, dimension(nkgrid,nee) :: iifun

end module global