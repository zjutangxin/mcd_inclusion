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

end module global