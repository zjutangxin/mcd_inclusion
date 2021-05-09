module routines
    implicit none
    integer, parameter, private :: dp = kind(1.0d0)

    contains

    function FFindInt(xx,xxVec,yyVec,nn)
        ! one-dimension interpolation
            implicit none
            integer, intent(in) :: nn
    
            real(dp), intent(in) :: xx
            real(dp), dimension(nn) :: xxvec, yyvec
            real(dp) :: FFindInt
    
            integer :: ind
            real(dp) :: slope
    
            if (xx <= xxVec(1)) then
                slope    = (yyVec(2)-yyVec(1))/(xxVec(2)-xxVec(1))
                FFindInt = yyVec(1)+slope*(xx-xxVec(1))
            elseif (xx >= xxVec(nn)) then
                slope    = (yyVec(nn)-yyVec(nn-1))/(xxVec(nn)-xxVec(nn-1))
                FFindInt = yyVec(nn)+slope*(xx-xxVec(nn))
            else
                ind = 1
                do while (xx > xxVec(ind))
                    ind = ind+1
                end do
                slope    = (yyVec(ind)-yyVec(ind-1))/(xxVec(ind)-xxVec(ind-1))
                FFindInt = yyVec(ind-1)+slope*(xx-xxVec(ind-1))
            endif
        end function FFindInt    

        subroutine FindIndex(xx,xxVec,xxGrid,ii,pp)

            implicit none
    
            integer, intent(in) :: xxGrid
            integer, intent(out) :: ii
            real(dp), intent(in) :: xx
            real(dp), dimension(xxGrid), intent(in) :: xxVec
            real(dp), intent(out) :: pp
    
            if (xx <= xxVec(1)) then
                ii = 2
                pp = 1.0_dp
            elseif (xx >= xxVec(xxGrid)) then
                ii = xxGrid
                pp = 0.0_dp
            else
                ii = 1
                do while (xx > xxVec(ii))
                    ii = ii+1
                end do
                pp = (xxVec(ii)-xx)/(xxVec(ii)-xxVec(ii-1))
            endif        
            
        end subroutine FindIndex        

        subroutine piksr2(n,arr,brr)

            implicit none
    
            integer, intent(in) :: n
            real(dp), dimension(n), intent(inout) :: arr(n), brr(n)
    
            integer :: i,j
            real(dp) :: a,b
    
            DO  J=2,N
                A=ARR(J)
                B=BRR(J)
                DO I=J-1,1,-1
                    IF (ARR(I) .LE. A) GO TO 10
                        ARR(I+1)=ARR(I)
                        BRR(I+1)=BRR(I)
                END DO
                I=0
                10    ARR(I+1)=A
                BRR(I+1)=B
            END DO        
        end subroutine piksr2
    
        function FindPerc(perc,a,p,n)
            ! Share of a by top perc in p
            implicit none
    
            integer, intent(in) :: n
            real(dp), dimension(n), intent(in) :: a,p
            real(dp), intent(in) :: perc
            real(dp) :: FindPerc
    
            real(dp), dimension(n) :: xcum,pcum
            integer :: indi
        
            pcum(1) = p(1)
            xcum(1) = a(1)*p(1)
    
            do indi = 2,n
                pcum(indi) = pcum(indi-1)+p(indi)
                xcum(indi) = xcum(indi-1)+a(indi)*p(indi)
            end do
            pcum = pcum/pcum(n)
            xcum = xcum/xcum(n)
    
            FindPerc = 1.0_dp - FFindInt(perc,pcum,xcum,n)
    
        end function FindPerc

end module routines