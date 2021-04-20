module parameters

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    real(dp), parameter :: beta = 0.9040_dp
    real(dp), parameter :: sigma = 1.5_dp
    real(dp), parameter :: alpha = 1.0_dp/3.0_dp
    real(dp), parameter :: delta = 0.06_dp
    real(dp), parameter :: nu = 0.21_dp
    real(dp), parameter :: eta = 4.15_dp
    real(dp), parameter :: psi = 0.8940_dp
    real(dp), parameter :: phi = 0.3670_dp

    real(dp), parameter :: lambda = 1.35_dp
    real(dp), parameter :: tplus = 0.50_dp
    real(dp), parameter :: tminus = -0.2975_dp
    real(dp), parameter :: q = 7.0_dp

    real(dp), parameter :: kmax = 12500_dp
    real(dp), parameter :: kmin = 1e-9
    real(dp), parameter :: kint = 25.0_dp

    integer, parameter :: nkgrid = 501
    integer, parameter :: sr = 10
    integer, parameter :: nkgridc = sr*(nkgrid-1)+1
    integer, parameter :: srsim = 15
    integer, parameter :: nsim = srsim*(nkgrid-1)+1
    integer, parameter :: nintk = 201
    integer, parameter :: negrid = 12
    integer, parameter :: nee = negrid*2
    
    real(dp), parameter :: tolv = 1e-6
    real(dp), parameter :: told = 1e-7
    real(dp), parameter :: tolmin = 1e-2
    integer, parameter :: maxiterv = 5000
    integer, parameter :: maxiterd = 10000
    
    ! minpack
    integer, parameter :: nmarket = 2
    integer, parameter :: lwa = (nmarket*(3*nmarket+13))/2

    ! MPI
    integer, parameter :: root = 0
    integer, parameter :: ni_total = nkgrid*nee


end module parameters