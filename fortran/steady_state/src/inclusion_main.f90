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
!
! Compiling Environment:
!   GNU gfortran on Ubuntu 16.04
!
! Library Used:
!   - N/A 
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
    integer :: iflag, info
    real(dp), dimension(lwa) :: wa
    external :: fcn_ss_cont

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

    write (*,*) 'Pareto distribution: compare mean and variance'
    write (*,*) 'Theoretical: ', emean, evar
    write (*,*) 'Numerical', nummean, numvar

    wguess = 0.0803852580114_dp
    rguess = -0.0501250_dp
    x_in(1) = wguess
    x_in(2) = rguess

    ! evaluate excess demand at price xin
    ! call fcn_ss(x_in,fvec)
    ! call fcn_ss_cont(nmarket,x_in,fvec,iflag)
    call hybrd1(fcn_ss_cont,nmarket,x_in,fvec,tolmin,info,wa,lwa)

    !write (*,*) 'wage == ', wguess, 'interest == ', rguess
    write (*,*) 'wage == ', x_in(1), 'interest == ', x_in(2)
    write (*,*) 'excess demand'
    write (*,*) 'labor == ', fvec(1), 'capital == ', fvec(2)

    call save_results
    call cpu_time(etime)

    write (*,*) 'time running == ', etime-btime, ' secs'

end program inclusion