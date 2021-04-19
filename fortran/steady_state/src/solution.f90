subroutine fcn_ss(x_in,fvec)

    use parameters
    use global
    implicit none

    real(dp), dimension(2), intent(in) :: x_in
    real(dp), dimension(2), intent(out) :: fvec
    real(dp), dimension(nkgrid,nee) :: vnn, income, vh
    real(dp), dimension(nsim,nee) :: distssn
    real(dp), dimension(nsim) :: disttemp
    real(dp), dimension(nkgridc) :: vtemp
    real(dp), dimension(negrid,2) :: lstr, kstr

    real(dp) :: w, r, kappa, mop, kcstr, lcstr, pistr, incomes
    real(dp) :: cons, vtmp
    real(dp) :: eps, epsh, tolh
    integer :: loopn, indz, inde, indk, indkc, maxiterh, looph
    integer :: jj, aind, nj, loopd, indj, induse

    w = x_in(1)
    r = x_in(2)

    write (*,*) 'r == ', r, ' w ==', w
    if (r+delta .le. 0.0_dp) then
        r = -delta + 1e-6
        write (*,*) 'warning!!!! negative kappa!!!!'
    end if

    kappa = alpha*w/((1-alpha)*(r+delta))
    lstr = (w/(tauz*(1-nu)*(1-alpha)*kappa**(alpha*(1-nu))))**(-1/nu)
    kstr = kappa*lstr

    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            kcstr = min(kstr(inde,indz),lambda*kgrid(indk))
            lcstr = (w/(tauz(inde,indz)*(1-nu)*(1-alpha)* &
                (kcstr)**(alpha*(1-nu))))**(1/(alpha*(nu-1)-nu))
            pistr = tauz(inde,indz)*(kcstr**alpha*lcstr**(1-alpha))**(1-nu) &
                - w*lcstr - (r+delta)*kcstr
            mop = max(pistr,w)

            if (pistr .ge. w ) then
                polk(indk,(indz-1)*negrid+inde) = kcstr
                poll(indk,(indz-1)*negrid+inde) = lcstr
                poly(indk,(indz-1)*negrid+inde) = &
                    egrid(inde)*(kcstr**alpha*lcstr**(1-alpha))**(1-nu)
                polpi(indk,(indz-1)*negrid+inde) = pistr
            end if
            income(indk,(indz-1)*negrid+inde) = mop + (1+r)*kgrid(indk)
            vn(indk,(indz-1)*negrid+inde) =(income(indk,(indz-1)*negrid+inde)**(1-sigma)-1)/(1-sigma)
        end do ! indk
    end do ! inde
    end do ! indz

    eps = 10000.0_dp
    loopn = 1

    do loopn = 1,maxiterv
    !if (mod(loopn,20) .eq. 0) then
        write (*,*) 'loopn == ', loopn, ' eps == ', eps
    !end if
    polaind = nkgridc 
    polc = 0.0_dp
    pola = 0.0_dp
    vtemp = 0.0_dp

    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            incomes = income(indk,(indz-1)*negrid+inde)
            if (incomes .le. kint) then
                nj = max(1,int((incomes-kmin)*sr*(nintk-1.0_dp)/(kint-kmin)+1.0_dp)-1)
            else
                nj = min(nkgridc,int((nintk-1)*sr+1+(incomes-kint)/(kmax-kint)*sr*(nkgrid-nintk)))-1    
            end if           
            nj = min(nj,polaind(min(nkgrid,indk+1),(indz-1)*negrid+inde))

            do indkc = nj,1,-1
                jj = min(nkgrid-1,(indkc-1)/sr+1)
                vtemp(indkc) = psi*&
                    ((vn(jj+1,inde+negrid*(indz-1))-vn(jj,inde+negrid*(indz-1))) &
                    /(kgrid(jj+1)-kgrid(jj))*&
                    (kgridc(indkc)-kgrid(jj))+vn(jj,inde+negrid*(indz-1)))
                vtemp(indkc) = vtemp(indkc) + (1-psi)*&
                    sum(omega*((vn(jj+1,:)-vn(jj,:))/(kgrid(jj+1)-kgrid(jj))*(kgridc(indkc)-kgrid(jj)) &
                    +vn(jj,:)))
                cons = incomes - kgridc(indkc)
                vtemp(indkc) = beta*vtemp(indkc) + (cons**(1-sigma)-1)/(1-sigma)

                if ((indkc .lt. nj) .and. (vtemp(indkc) .lt. vtemp(min(indkc+1,nj))) ) then
                    exit
                end if
                vnn(indk,inde+(indz-1)*negrid) = vtemp(indkc)
                pola(indk,inde+(indz-1)*negrid) = kgridc(indkc)
                polc(indk,inde+(indz-1)*negrid) = incomes - kgridc(indkc)
                polaind(indk,inde+(indz-1)*negrid) = indkc   

            end do ! indkc

        end do ! indk
    end do ! inde
    end do ! indz

    ! Howard acceleration
    maxiterh = maxiterv
    looph = 1
    epsh = 1
    vh = vnn
    tolh = tolv

    do looph = 1,maxiterv
    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            aind = polaind(indk,inde+(indz-1)*negrid)
            jj = min(nkgrid-1,int((aind-1)/sr)+1)
            vtmp = psi*&
                ((vh(jj+1,inde+negrid*(indz-1))-vh(jj,inde+negrid*(indz-1)))&
                /(kgrid(jj+1)-kgrid(jj))*&
                (kgridc(aind)-kgrid(jj))+vh(jj,inde+negrid*(indz-1)))
            vtmp = vtmp + (1-psi)*&
                sum(omega*((vh(jj+1,:)-vh(jj,:))/(kgrid(jj+1)-kgrid(jj))*(kgridc(aind)-kgrid(jj))&
                +vh(jj,:)))
            vnn(indk,inde+negrid*(indz-1)) = (polc(indk,inde+negrid*(indz-1))**(1-sigma)-1)/(1-sigma)&
                + beta*vtmp
        end do ! indk
    end do ! inde
    end do ! indz

    epsh = maxval(dabs(vnn - vh))
    vh = vnn 
    if (eps .lt. tolh) exit

    end do ! looph, howard

    eps = maxval(dabs(vn-vnn))
    vn = vnn 
    if (eps .lt. tolv) exit

    enddo ! loopn

    eps = 10000.0_dp
    dpolc = 0.0_dp
    dpola = 0.0_dp
    dpolk = 0.0_dp
    dpoll = 0.0_dp
    dpoly = 0.0_dp

    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nsim
            jj = min(nkgrid-1,int((indk-1)/srsim)+1) 
            dpolc(indk,inde+(indz-1)*negrid) = &
                ((polc(jj+1,inde+negrid*(indz-1))-polc(jj,inde+negrid*(indz-1)))&
                /(kgrid(jj+1)-kgrid(jj))*&
                (kgrids(indk)-kgrid(jj))+polc(jj,inde+negrid*(indz-1)));        
            dpola(indk,inde+(indz-1)*negrid) = &
                ((pola(jj+1,inde+negrid*(indz-1))-pola(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+pola(jj,inde+negrid*(indz-1)))
            dpolk(indk,inde+(indz-1)*negrid) = &
                 ((polk(jj+1,inde+negrid*(indz-1))-polk(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+polk(jj,inde+negrid*(indz-1)))
            dpoll(indk,inde+(indz-1)*negrid) = &
                ((poll(jj+1,inde+negrid*(indz-1))-poll(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+poll(jj,inde+negrid*(indz-1)))
            dpoly(indk,inde+(indz-1)*negrid) = &
                ((poly(jj+1,inde+negrid*(indz-1))-poly(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+poly(jj,inde+negrid*(indz-1)))
            dpolpi(indk,inde+(indz-1)*negrid) = &
                ((polpi(jj+1,inde+negrid*(indz-1))-poly(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+polpi(jj,inde+negrid*(indz-1)))            
        end do
    end do
    end do

    ddistss = 1.0_dp/nsim
    do loopd = 1,maxiterd
        if (mod(loopd,20) .eq. 0) then
            write (*,*) 'loopd == ', loopd, ' eps == ', eps
        end if

        distssn = 0.0_dp
        do indz = 1,2
        do inde = 1,negrid
            do indk = 1,nsim
                if (dpola(indk,inde+(indz-1)*negrid) .le. kint) then
                    jj = max(1,int((dpola(indk,inde+(indz-1)*negrid)-kmin)*srsim*(nintk-1.0)/(kint-kmin)+1.0))
                else
                    jj = min(nsim-1,int((nintk-1)*srsim+1+&
                        (dpola(indk,inde+(indz-1)*negrid)-kint)/(kmax-kint)*srsim*(nkgrid-nintk)))
                end if
                distssn(jj,inde+(indz-1)*negrid) = distssn(jj,inde+(indz-1)*negrid) + &
                    ddistss(indk,inde+(indz-1)*negrid)*&
                    (kgrids(jj+1)-dpola(indk,inde+(indz-1)*negrid))/(kgrids(jj+1)-kgrids(jj))
                distssn(jj+1,inde+(indz-1)*negrid) = distssn(jj+1,inde+(indz-1)*negrid) + &
                    ddistss(indk,inde+(indz-1)*negrid)*&
                    (dpola(indk,inde+(indz-1)*negrid)-kgrids(jj))/(kgrids(jj+1)-kgrids(jj))
            enddo ! indk
        enddo ! inde
        enddo ! indz

        disttemp = matmul(ddistss,omega)

        do indj = 1,nee
            distssn(:,indj) = psi*distssn(:,indj) + (1-psi)*disttemp
        end do

        eps = maxval(dabs(distssn-ddistss))
        ddistss = distssn
        if (eps .lt. told) exit
    enddo

    ! open(1,file='./results/distend.txt',form='formatted')
    ! write(1,'(99f16.6)') ddistss(nsim,:)
    ! close(1)
    
    ! compute aggregate demand and supply
    asupply = sum(matmul(dpola*ddistss,omega))
    ademand = sum(matmul(dpolk*ddistss,omega))
    ldemand = sum(matmul(dpoll*ddistss,omega))
    lsupply = 0.0_dp
    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nsim
            if (dpoll(indk,inde+(indz-1)*negrid) .eq. 0.0_dp) then
                lsupply = lsupply + &
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
            end if
        enddo ! indk
    enddo ! inde
    enddo ! indz

    fvec(1) = ldemand - lsupply
    fvec(2) = ademand - asupply

end subroutine fcn_ss


subroutine fcn_ss_cont(np,x_in,fvec,iflag)

    ! aggregate labor supply is interpolated

    use parameters
    use global
    implicit none

    integer, intent(in) :: np
    integer, intent(out) :: iflag
    real(dp), dimension(2), intent(in) :: x_in
    real(dp), dimension(2), intent(out) :: fvec
    real(dp), dimension(nkgrid,nee) :: vnn, income, vh
    real(dp), dimension(nee) :: wgt, maxprof
    real(dp), dimension(nsim,nee) :: distssn
    real(dp), dimension(nsim) :: disttemp
    real(dp), dimension(nkgridc) :: vtemp
    real(dp), dimension(negrid,2) :: lstr, kstr

    real(dp) :: w, r, kappa, mop, kcstr, lcstr, pistr, incomes
    real(dp) :: cons, vtmp
    real(dp) :: eps, epsh, tolh
    integer :: loopn, indz, inde, indk, indkc, maxiterh, looph
    integer :: jj, aind, nj, loopd, indj, induse

    w = x_in(1)
    r = x_in(2)

    write (*,*) 'r == ', r, ' w ==', w
    if (r+delta .le. 0.0_dp) then
        r = -delta + 1e-6
        write (*,*) 'warning!!!! negative kappa!!!!'
    end if

    kappa = alpha*w/((1-alpha)*(r+delta))
    lstr = (w/(tauz*(1-nu)*(1-alpha)*kappa**(alpha*(1-nu))))**(-1/nu)
    kstr = kappa*lstr

    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            kcstr = min(kstr(inde,indz),lambda*kgrid(indk))
            lcstr = (w/(tauz(inde,indz)*(1-nu)*(1-alpha)* &
                (kcstr)**(alpha*(1-nu))))**(1/(alpha*(nu-1)-nu))
            pistr = tauz(inde,indz)*(kcstr**alpha*lcstr**(1-alpha))**(1-nu) &
                - w*lcstr - (r+delta)*kcstr
            mop = max(pistr,w)

            polk(indk,(indz-1)*negrid+inde) = kcstr
            poll(indk,(indz-1)*negrid+inde) = lcstr
            poly(indk,(indz-1)*negrid+inde) = &
                egrid(inde)*(kcstr**alpha*lcstr**(1-alpha))**(1-nu)
            polpi(indk,(indz-1)*negrid+inde) = pistr

            if (pistr .ge. w ) then
                polo(indk,(indz-1)*negrid+inde) = 1 ;
            end if
            income(indk,(indz-1)*negrid+inde) = mop + (1+r)*kgrid(indk)
            vn(indk,(indz-1)*negrid+inde) =(income(indk,(indz-1)*negrid+inde)**(1-sigma)-1)/(1-sigma)
        end do ! indk
    end do ! inde
    end do ! indz

    eps = 10000.0_dp
    loopn = 1

    do loopn = 1,maxiterv
    if (mod(loopn,5) .eq. 0) then
        write (*,*) 'loopn == ', loopn, ' eps == ', eps
    end if
    polaind = nkgridc 
    polc = 0.0_dp
    pola = 0.0_dp
    vtemp = 0.0_dp

    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            incomes = income(indk,(indz-1)*negrid+inde)
            if (incomes .le. kint) then
                nj = max(1,int((incomes-kmin)*sr*(nintk-1.0_dp)/(kint-kmin)+1.0_dp)-1)
            else
                nj = min(nkgridc,int((nintk-1)*sr+1+(incomes-kint)/(kmax-kint)*sr*(nkgrid-nintk)))-1    
            end if           
            nj = min(nj,polaind(min(nkgrid,indk+1),(indz-1)*negrid+inde))

            do indkc = nj,1,-1
                jj = min(nkgrid-1,(indkc-1)/sr+1)
                vtemp(indkc) = psi*&
                    ((vn(jj+1,inde+negrid*(indz-1))-vn(jj,inde+negrid*(indz-1))) &
                    /(kgrid(jj+1)-kgrid(jj))*&
                    (kgridc(indkc)-kgrid(jj))+vn(jj,inde+negrid*(indz-1)))
                vtemp(indkc) = vtemp(indkc) + (1-psi)*&
                    sum(omega*((vn(jj+1,:)-vn(jj,:))/(kgrid(jj+1)-kgrid(jj))*(kgridc(indkc)-kgrid(jj)) &
                    +vn(jj,:)))
                cons = incomes - kgridc(indkc)
                vtemp(indkc) = beta*vtemp(indkc) + (cons**(1-sigma)-1)/(1-sigma)

                if ((indkc .lt. nj) .and. (vtemp(indkc) .lt. vtemp(min(indkc+1,nj))) ) then
                    exit
                end if
                vnn(indk,inde+(indz-1)*negrid) = vtemp(indkc)
                pola(indk,inde+(indz-1)*negrid) = kgridc(indkc)
                polc(indk,inde+(indz-1)*negrid) = incomes - kgridc(indkc)
                polaind(indk,inde+(indz-1)*negrid) = indkc   

            end do ! indkc

        end do ! indk
    end do ! inde
    end do ! indz

    ! Howard acceleration
    maxiterh = maxiterv
    looph = 1
    epsh = 1
    vh = vnn
    tolh = tolv

    do looph = 1,maxiterv
    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nkgrid
            aind = polaind(indk,inde+(indz-1)*negrid)
            jj = min(nkgrid-1,int((aind-1)/sr)+1)
            vtmp = psi*&
                ((vh(jj+1,inde+negrid*(indz-1))-vh(jj,inde+negrid*(indz-1)))&
                /(kgrid(jj+1)-kgrid(jj))*&
                (kgridc(aind)-kgrid(jj))+vh(jj,inde+negrid*(indz-1)))
            vtmp = vtmp + (1-psi)*&
                sum(omega*((vh(jj+1,:)-vh(jj,:))/(kgrid(jj+1)-kgrid(jj))*(kgridc(aind)-kgrid(jj))&
                +vh(jj,:)))
            vnn(indk,inde+negrid*(indz-1)) = (polc(indk,inde+negrid*(indz-1))**(1-sigma)-1)/(1-sigma)&
                + beta*vtmp
        end do ! indk
    end do ! inde
    end do ! indz

    epsh = maxval(dabs(vnn - vh))
    vh = vnn 
    if (eps .lt. tolh) exit

    end do ! looph, howard

    eps = maxval(dabs(vn-vnn))
    vn = vnn 
    if (eps .lt. tolv) exit

    enddo ! loopn

    eps = 10000.0_dp
    dpolc = 0.0_dp
    dpola = 0.0_dp
    dpolk = 0.0_dp
    dpoll = 0.0_dp
    dpoly = 0.0_dp

    do indz = 1,2
    do inde = 1,negrid
        do indk = 1,nsim
            jj = min(nkgrid-1,int((indk-1)/srsim)+1) 
            dpolc(indk,inde+(indz-1)*negrid) = &
                ((polc(jj+1,inde+negrid*(indz-1))-polc(jj,inde+negrid*(indz-1)))&
                /(kgrid(jj+1)-kgrid(jj))*&
                (kgrids(indk)-kgrid(jj))+polc(jj,inde+negrid*(indz-1)));        
            dpola(indk,inde+(indz-1)*negrid) = &
                ((pola(jj+1,inde+negrid*(indz-1))-pola(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+pola(jj,inde+negrid*(indz-1)))
            dpolk(indk,inde+(indz-1)*negrid) = &
                 ((polk(jj+1,inde+negrid*(indz-1))-polk(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+polk(jj,inde+negrid*(indz-1)))
            dpoll(indk,inde+(indz-1)*negrid) = &
                ((poll(jj+1,inde+negrid*(indz-1))-poll(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+poll(jj,inde+negrid*(indz-1)))
            dpoly(indk,inde+(indz-1)*negrid) = &
                ((poly(jj+1,inde+negrid*(indz-1))-poly(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+poly(jj,inde+negrid*(indz-1)))
            dpolpi(indk,inde+(indz-1)*negrid) = &
                ((polpi(jj+1,inde+negrid*(indz-1))-polpi(jj,inde+negrid*(indz-1)))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+polpi(jj,inde+negrid*(indz-1)))            
            dpolo(indk,inde+(indz-1)*negrid) = &
                nint( (real(polo(jj+1,inde+negrid*(indz-1)))-real(polo(jj,inde+negrid*(indz-1))))&
                        /(kgrid(jj+1)-kgrid(jj))*&
                        (kgrids(indk)-kgrid(jj))+real(polo(jj,inde+negrid*(indz-1))) )
        end do
    end do
    end do

    ddistss = 1.0_dp/nsim
    do loopd = 1,maxiterd
        if (mod(loopd,50) .eq. 0) then
            write (*,*) 'loopd == ', loopd, ' eps == ', eps
        end if

        distssn = 0.0_dp
        do indz = 1,2
        do inde = 1,negrid
            do indk = 1,nsim
                if (dpola(indk,inde+(indz-1)*negrid) .le. kint) then
                    jj = max(1,int((dpola(indk,inde+(indz-1)*negrid)-kmin)*srsim*(nintk-1.0)/(kint-kmin)+1.0))
                else
                    jj = min(nsim-1,int((nintk-1)*srsim+1+&
                        (dpola(indk,inde+(indz-1)*negrid)-kint)/(kmax-kint)*srsim*(nkgrid-nintk)))
                end if
                distssn(jj,inde+(indz-1)*negrid) = distssn(jj,inde+(indz-1)*negrid) + &
                    ddistss(indk,inde+(indz-1)*negrid)*&
                    (kgrids(jj+1)-dpola(indk,inde+(indz-1)*negrid))/(kgrids(jj+1)-kgrids(jj))
                distssn(jj+1,inde+(indz-1)*negrid) = distssn(jj+1,inde+(indz-1)*negrid) + &
                    ddistss(indk,inde+(indz-1)*negrid)*&
                    (dpola(indk,inde+(indz-1)*negrid)-kgrids(jj))/(kgrids(jj+1)-kgrids(jj))
            enddo ! indk
        enddo ! inde
        enddo ! indz

        disttemp = matmul(ddistss,omega)

        do indj = 1,nee
            distssn(:,indj) = psi*distssn(:,indj) + (1-psi)*disttemp
        end do

        eps = maxval(dabs(distssn-ddistss))
        ddistss = distssn
        if (eps .lt. told) exit
    enddo
   
    ! compute aggregate demand and supply
    ! asupply = sum(matmul(dpola*ddistss,omega))
    ! ademand = sum(matmul(dpolk*ddistss,omega))
    ! ldemand = sum(matmul(dpoll*ddistss,omega))
    ! lsupply = 0.0_dp
    ! do indz = 1,2
    ! do inde = 1,negrid
    !     do indk = 1,nsim
    !         if (dpoll(indk,inde+(indz-1)*negrid) .eq. 0.0_dp) then
    !             lsupply = lsupply + &
    !                 ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
    !         end if
    !     enddo ! indk
    ! enddo ! inde
    ! enddo ! indz

    ! compute aggregate demand and supply
    ! with aggregate labor supply continuous
    asupply = 0.0_dp
    ademand = 0.0_dp
    lsupply = 0.0_dp
    ldemand = 0.0_dp
    wgt = 0.0_dp 
    maxprof = maxval(polpi,1)

    if (maxprof(1) > w) then
        wgt(1) = w/maxprof(1)
    end if
    if (maxprof(negrid+1) > w) then
        wgt(negrid+1) = w/maxprof(negrid+1)
    end if

    do indz = 1,2
    do inde = 1,negrid-1
        induse = inde+(indz-1)*negrid
        if ((w > maxprof(induse)) .and. (w > maxprof(induse+1))) then
            wgt(induse+1) = 1.0_dp
        else if ((w > maxprof(induse)) .and. (w < maxprof(induse+1))) then
            wgt(induse+1) = (w-maxprof(induse))/(maxprof(induse+1)-maxprof(induse))
        end if
    end do
    end do

    do indz = 1,2
    do inde = 1,negrid
        if ((wgt(inde+(indz-1)*negrid)>0) .and. (wgt(inde+(indz-1)*negrid) < 1)) then
            do indk = 1,nsim
                if (dpolo(indk,inde+(indz-1)*negrid) .gt. 0) then
                    ademand = ademand + dpolk(indk,inde+(indz-1)*negrid)*(1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    ldemand = ldemand + dpoll(indk,inde+(indz-1)*negrid)*(1-wgt(inde+(indz-1)*negrid))*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    lsupply = lsupply + ddistss(indk,inde+(indz-1)*negrid)&
                        *omega(inde+(indz-1)*negrid)*wgt(inde+(indz-1)*negrid)
                else
                    asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    lsupply = lsupply + &
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                end if
            end do
        else
            do indk = 1,nsim
                if (dpolo(indk,inde+(indz-1)*negrid) .gt. 0) then
                    ademand = ademand + dpolk(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    ldemand = ldemand + dpoll(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                else
                    asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*&
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                    lsupply = lsupply + &
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)
                end if
            end do
        end if
    end do
    end do

    fvec(1) = ldemand - lsupply
    fvec(2) = ademand - asupply
end subroutine fcn_ss_cont
