module ModuleSaving
     
use globals
    
    contains
    
subroutine SaveSteadys    
    implicit none
    
    open(1,  file = 'Output/phi.txt',      status = 'unknown')
    open(2,  file = 'Output/aplus_pol.txt',  status = 'unknown')
    open(3,  file = 'Output/c.txt',  status = 'unknown')
    open(4,  file = 'Output/V.txt',  status = 'unknown')        
    open(5,  file = 'Output/RHS.txt',      status = 'unknown')
    open(6,  file = 'Output/pi.txt',  status = 'unknown')
 

    
    write(1, *) phi
    write(2, *) aplus_pol
    write(3, *) c
    write(4, *) V     
    write(5,*) RHS
    write(6,*) pi

    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
  
end subroutine SaveSteadys

end module  ModuleSaving
