subroutine save_results

    use global
    use parameters
    implicit none

    integer :: indi

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

    open(1,file='./results/e.txt',form='formatted')
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

end subroutine save_results