
    function FOCnewton_wa(tla,tlw,tlv) result(c_out)
    
    implicit none
    
    real(8), intent(in)             :: tla(:), tlw(:), tlv
    real(8), dimension(size(tla,1)) :: caux, cnew, fval, dfval, c_out
    
    real(8) :: tol_dif, err_dif
    integer :: itfoc, maxitfoc, itbnd 
    
    tol_dif = 1e-9; maxitfoc = 250; err_dif = 1.0D0
    
    caux = tiny
    
    itfoc = 1
    
    do while ( err_dif > tol_dif .and. itfoc <= maxitfoc )
        
        fval  = tlw*(caux**(-tlv)) - ( caux + tla )
        
        dfval = - (tlw*tlv*(caux**(-tlv-1.0D0))) - 1.0D0 
        
        err_dif = maxval( abs(fval) )
        
        where (abs(fval) > tol_dif) caux = caux - (fval/dfval)        
        
        itfoc = itfoc + 1
        
        ! print *, 'err_dif = ', err_dif, 'itfoc = ', itfoc 
        
        if (itfoc == maxitfoc - 1) then
            print *, 'couldnt solve the bisection'
            print *, 'err_dif = ', err_dif
        end if        
    end do    
    
    c_out = caux
    
    end function FOCnewton_wa