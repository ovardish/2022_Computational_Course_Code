!##############################################################################
! Oliko Vardishvili
! oliko.vardishvili@yale.edu
!
!##############################################################################
module globals


    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.50d0
     real*8, parameter :: egam = 1d0 - 1d0/gamma 
    real*8, parameter :: beta = 0.96d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 0.08d0
    real*8, parameter :: rho = 0.2d0
    real*8, parameter :: sigma_eps = 0.4d0*sqrt(1d0-rho**2)

    ! numerical parameters
    real*8, parameter :: a_l = 0.0d0
    real*8, parameter :: a_u = 25.0d0
    real*8, parameter :: a_grow = 0.01d0
    real*8, parameter :: sig_in = 1d-10
    real*8, parameter :: sig_out = 1d-6

    
    integer, parameter :: itermax = 50000

    ! macroeconomic variables
    real*8 :: r, w, KK, AA, LL, YY, CC, II

    ! the shock process
    integer, parameter :: NS_p = 7
    integer, parameter :: NS = 7  ! NS_p+1
    real*8 :: pi_p(NS_p, NS_p), eta_p(NS_p), weights_p(NS_p)
    real*8 :: pi(NS, NS), eta(NS), weights(NS)

    ! policy and value function
    integer, parameter :: NA = 1002
    real*8 :: a(0:NA), c_pol(0:NA, NS), V(0:NA, NS)

    ! variables to numerically determine policy function
    real*8  ::  c_new(0:NA, NS), V_new(0:NA, NS), RHS(0:NA, NS) , EV(0:NA, NS)
    real*8  ::  a_pol(0:NA, NS)
    logical ::  check

    ! variables to communicate with function
    real*8  :: a_com
    integer :: is_com
    integer :: ia_com



    ! variables to numerically determine the steady state distribution
    real*8 :: phi(0:NA, NS), phi_new(0:NA, NS)
    
    ! Value function interpolation
    
    real*8 :: coeff_V(NA+3, NS), coeff_c(NA+3, NS)
    real*8 :: con_lev, x_in, fret

    integer :: maximization_iter
    integer :: grid_method
    contains



          
  
                   

    
    
 
    
    
    
    
    
    
    
    
    
    
    
    
    
          ! the function that should be minimized
    function utility(x_in)

        use toolbox
        
        implicit none
        real*8, intent(in) :: x_in
        real*8 :: utility, cons, vplus , varphi
        integer :: is_p, ial, iar
        
       ! maximization_iter=maximization_iter+1    
       ! print*, 'maximization_iter=', maximization_iter
        !pause
        
        ! calculate consumption
        cons = max(a(ia_com)*(1+r) + eta(is_com)*w- x_in, 1d-10)
        
        call linint_Grow(x_in, a_l, a_u, a_grow, NA, ial, iar, varphi)
        ! get utility function
        utility = - (cons**egam/egam + beta*max(varphi*EV(ial,  is_com) + (1d0-varphi)*EV( iar, is_com), 1d-10)**egam/egam)

    end function  
     
    
       ! calculates the invariant distribution of households
    subroutine get_distribution()
        use toolbox
        implicit none
        integer :: ia, is, iter, is_p
        real*8 :: aplus, con_lev
        integer :: ial(0:NA, NS), iar(0:NA, NS)
        real*8 :: varphi(0:NA, NS)

         
        ! get the interpolation shares and points
        do ia = 0, NA
            do is = 1, NS

                ! calculate where this guy would go
               aplus = (1d0+r)*a(ia) + w*eta(is) - c_pol(ia,is)
               aplus = min(max(aplus, a_l), a_u)
                ! determine the gridpoints in between this decision lies
              
              call linint_grow(aplus, a_l, a_u, a_grow,  NA, ial(ia, is), iar(ia, is), varphi(ia, is))
            enddo
        enddo
        

        ! iterate until the distribution function converges
        do iter = 1, itermax

            phi_new = 0d0

            do ia = 0, NA
                do is = 1, NS

                    do is_p = 1, NS
                        phi_new(ial(ia,is), is_p) = phi_new(ial(ia,is), is_p) + &
                            pi(is,is_p)*varphi(ia,is)*phi(ia,is)
                        phi_new(iar(ia,is), is_p) = phi_new(iar(ia,is), is_p) + &
                            pi(is,is_p)*(1d0-varphi(ia,is))*phi(ia,is)
                    enddo
                enddo
            enddo

            con_lev = maxval(abs(phi_new(:, :) - phi(:, :))/max(abs(phi(:, :)), 1d-10))

            ! update distribution
            phi = phi_new

            ! check for convergence
            if(con_lev < sig_in)then
                 if (abs(sum(phi)-1.0d0) > 0.000001d0)    print*, 'your distribution goes nuts'
                return
            endif
        enddo

        write(*,*)'Distribution Function Iteration did not converge'

    end subroutine
    

    
    
    
    
    
    
    

                   
end module