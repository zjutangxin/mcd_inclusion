subroutine compute_moments

    use parameters
    use global
    use routines
    implicit none
    integer :: indz, inde, indk, nldvec, nidvec, index
    real(dp) :: perc
    real(dp), dimension(:), allocatable :: ld_vec, id_vec, dld_vec, did_vec

    ! --------------- Consistency check ----------------
    ! government budget
    govtax = 0.0_dp
    govsub = 0.0_dp
    do inde = 1,negrid
        if ((wgt(inde) .gt. 0.0_dp) .and. (wgt(inde) .lt. 1.0_dp)) then
            do indk = 1,nsim
                if (dpolo(indk,inde) .gt. 0) then
                    govtax = govtax + dpoly(indk,inde)*(1-wgt(inde))*&
                        ddistss(indk,inde)*omega(inde)*tplus
                end if
            end do
        else
            do indk = 1,nsim
                if (dpolo(indk,inde) .gt. 0) then
                    govtax = govtax + dpoly(indk,inde)*&
                        ddistss(indk,inde)*omega(inde)*tplus
                end if
            end do
        end if
    end do

    do inde = 1,negrid
        if ((wgt(inde+negrid).gt. 0.0_dp) .and. (wgt(inde+negrid) .lt. 1.0_dp)) then
            do indk = 1,nsim
                if (dpolo(indk,inde+negrid) .gt. 0) then
                    govsub = govsub + dpoly(indk,inde+negrid)*(1-wgt(inde+negrid))*&
                        ddistss(indk,inde+negrid)*omega(inde+negrid)*tminus
                end if
            end do
        else
            do indk = 1,nsim
                if (dpolo(indk,inde+negrid) .gt. 0) then
                    govsub = govsub + dpoly(indk,inde+negrid)*&
                        ddistss(indk,inde+negrid)*omega(inde+negrid)*tminus
                end if
            end do
        end if
    end do

    ! macro aggregates
    aggy = 0.0_dp
    aggc = 0.0_dp
    aggk = ademand
    aggef = 0.0_dp
    aggl = ldemand
    nfirm = 0.0_dp
    emeanf = 0.0_dp
    pimean = 0.0_dp
    incmean = 0.0_dp

    do indz = 1,2
    do inde = 1,negrid
        if ((wgt(inde+(indz-1)*negrid).gt. 0.0_dp) .and. (wgt(inde+(indz-1)*negrid) .lt. 1.0_dp)) then
            do indk = 1,nsim
                if ( dpolo(indk,inde+(indz-1)*negrid) .gt. 0 ) then
                    aggy = aggy + dpoly(indk,inde+(indz-1)*negrid)*(1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    aggc = aggc + dpolc(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    aggef = aggef + (dpolk(indk,inde+(indz-1)*negrid)-kgrids(indk))*&
                        (1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    nfirm = nfirm + (1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    emeanf = emeanf + egrid(inde)*(1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    pimean = pimean + dpolpi(indk,inde+(indz-1)*negrid)*(1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                else
                    aggc = aggc + dpolc(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                end if
            end do
        else
            do indk = 1,nsim
                if ( dpolo(indk,inde+(indz-1)*negrid) .gt. 0 ) then
                    aggy = aggy + dpoly(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    aggc = aggc + dpolc(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    aggef = aggef + (dpolk(indk,inde+(indz-1)*negrid)-kgrids(indk))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    nfirm = nfirm + ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    emeanf = emeanf + egrid(inde)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    pimean = pimean + dpolpi(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                else
                    aggc = aggc + dpolc(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                end if                
            end do
        end if
    end do
    end do 

    aggef = aggef/aggy          ! external finance ratio
    lmean = aggl/nfirm          ! labor demand per firm
    yperfirm = aggy/nfirm       ! output per firm
    lprod = aggy/aggl           ! labor productivity
    pimean = pimean/nfirm       ! average profit
    incmean = pimean*(1-lsupply)+lsupply*wss ! average income
    emeanf = emeanf/nfirm       ! average productivity

    ! resource constraint
    govex = govtax + govsub
    rcex = aggy - aggc - aggk*delta - govex

    ! --------------- Calibration ----------------
    nldvec = sum(sum(dpolo,1))
    nidvec = nsim*nee
    allocate(ld_vec(nldvec), id_vec(nidvec),dld_vec(nldvec),did_vec(nidvec))
    
    index = 0 
    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nsim
            index = index + 1
            if (dpolo(indk,inde+(indz-1)*negrid) .gt. 0) then
                id_vec(index) = dpolpi(indk,inde+(indz-1)*negrid)
            else
                id_vec(index) = wss
            end if
            did_vec(index) = ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
        end do
    end do
    end do
    did_vec = did_vec/sum(did_vec)  

    index = 0
    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nsim
            if (dpolo(indk,inde+(indz-1)*negrid) .gt. 0) then 
                index = index + 1
                ld_vec(index) = dpoll(indk,inde+(indz-1)*negrid)
                dld_vec(index) = ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
            end if
        end do
    end do
    end do
    dld_vec = dld_vec/sum(dld_vec)
    
    ! Top 10% employment
    call piksr2(nldvec,ld_vec,dld_vec)
    perc = 0.90_dp
    share_ld = FindPerc(perc,ld_vec,dld_vec,nldvec)

    ! Top 5% earning
    call piksr2(nidvec,id_vec,did_vec)
    perc = 0.95_dp
    share_id = FindPerc(perc,id_vec,did_vec,nidvec)

    ! Firm Exit Rate
    fexit = (1-psi)*lsupply

    ! --------------- Macro moments ----------------
    ! Print screen
    write (*,*) '========================================================'
    write (*,*) 'Calibration Moments'
    write (*,*) '--------------------------------------------------------'
    write (*,106) 'Parameters', 'Values', 'Moments', 'Data', 'Model'
    write (*,107) 'beta', beta, 'int rate', 0.045, rss
    write (*,107) 'eta', eta, 'top 10% ld', 0.67, share_ld
    write (*,107) 'nu', nu, 'top 5% inc', 0.30, share_id
    write (*,107) 'psi', psi, 'exit rate', 0.1, fexit
    write (*,107) 'lambda', lambda, 'external fin', 0.6, aggef
    write (*,*) '========================================================'

    write (*,*) ''
    write (*,*) '========================================================'
    write (*,*) 'Consistency Check'
    write (*,*) '--------------------------------------------------------'
    write (*,108) 'govex == ', govex
    write (*,108) 'rcex == ', rcex
    write (*,*) '========================================================'

106 format(a12,a10,a15,a10,a10)
107 format(a12, f10.4, a15, 2f10.4)
108 format(a12, f12.6)
    deallocate(ld_vec, id_vec, dld_vec, did_vec)    

end subroutine compute_moments

subroutine save_results

    use global
    use parameters
    implicit none

    integer :: indi

    ! grid points and policy functions
    open(1,file='./results/kgrid.txt',form='formatted')
    do indi = 1,nkgrid
        write(1,'(20ES14.6)') kgrid(indi)
    end do
    101 format (f16.6)
    close(1)    

    open(1,file='./results/kgridc.txt',form='formatted')
    do indi = 1,nkgridc
        write(1,'(20ES14.6)') kgridc(indi)
    end do
    close(1)    

    open(1,file='./results/kgrids.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(20ES14.6)') kgrids(indi)
    end do
    close(1)    

    open(1,file='./results/egrid.txt',form='formatted')
    do indi = 1,negrid
        write(1,'(20ES14.6)') egrid(indi), edist(indi)
    end do
    close(1)        

    open(1,file='./results/vn.txt',form='formatted')
    do indi = 1,nkgrid
        write(1,'(99ES14.6)') vn(indi,:)
    end do
    close(1)        

    open(1,file='./results/dpola.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(99ES14.6)') dpola(indi,:)
    end do
    close(1)        

    open(1,file='./results/dpoll.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(99ES14.6)') dpoll(indi,:)
    end do
    close(1)        
    
    open(1,file='./results/dpolk.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(99ES14.6)') dpolk(indi,:)
    end do
    close(1)            
    
    open(1,file='./results/dpoly.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(99ES14.6)') dpoly(indi,:)
    end do
    close(1)            
    
    open(1,file='./results/dpolpi.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(99ES14.6)') dpolpi(indi,:)
    end do
    close(1)                
    
    open(1,file='./results/ddistss.txt',form='formatted')
    do indi = 1,nsim
        write(1,'(99F16.10)') ddistss(indi,:)
    end do
    close(1)            

    ! model parameters
    open(1,file='./results/param_name.txt',form='formatted')
        write (1,200) 'beta', &
                      'sigma', &
                      'alpha', &
                      'delta', &
                      'nu', &
                      'eta', &
                      'psi', &
                      'lambda', &
                      'tplus', &
                      'tminus', &
                      'q', &
                      'zscaler'
    200 format(100(a15 /))                      
    close(1)

    open(1,file='./results/param_val.txt',form='formatted')
        write (1,201) beta
        write (1,201) sigma
        write (1,201) alpha
        write (1,201) delta
        write (1,201) nu
        write (1,201) eta
        write (1,201) psi
        write (1,201) lambda
        write (1,201) tplus
        write (1,201) tminus
        write (1,201) q
        write (1,201) zscaler
    201 format(f16.10)                      
    close(1)

    ! macro moments
    open(1,file='./results/moments_name.txt',form='formatted')
        write (1,200) 'aggy', &
                      'aggc', &
                      'aggk', &
                      'aggef', &
                      'aggl', &
                      'nfirm', &
                      'emeanf', &
                      'pimean', &
                      'incmean', &
                      'lmean', &
                      'yperfirm', &
                      'lprod', &
                      'share_ld', &
                      'share_id', &
                      'fexit'
    close(1)

    open(1,file='./results/moments_val.txt',form='formatted')
        write (1,201) aggy
        write (1,201) aggc
        write (1,201) aggk
        write (1,201) aggef
        write (1,201) aggl
        write (1,201) nfirm
        write (1,201) emeanf
        write (1,201) pimean
        write (1,201) incmean
        write (1,201) lmean
        write (1,201) yperfirm
        write (1,201) lprod
        write (1,201) share_ld
        write (1,201) share_id
        write (1,201) fexit
    close(1)

end subroutine save_results