

program MAIN

    use globals
    use ModuleSaving
    
    implicit none
    real*8 :: guess
    integer ::it, ia, is, is_p
    
    !initialize variables
    call initialize()

    !set initial guess of the interest rate
    r=0.04                                     ! Initialize the interest rate
 
    
    ! solve_steady
    call solve_Steady()      
    
    ! Store policies in the benchmark economy
    call value_function                        ! calculate value function given policy functions
    
    

contains

 

    
    ! For solving the model
    subroutine solve_Steady

        use toolbox

        implicit none
        real*8 :: x
        logical :: check

        ! start timer
        call tic()

        write(*, '(a)')'       K/Y         r      diff_K'

        ! set tolerance level for outer optimization
        call settol_root(sig_out)

        ! set initial guess
        x = r
        check = .false.
     
        ! solve asset market equilibrium
        call fzero(x, asset_market, check)

        if(check)write(*,'(a)')'ERROR IN ROOTFINDING'

        call toc
     
        call SaveSteadys

    end subroutine


    ! For initializing variables
    subroutine initialize()

        use toolbox

        implicit none
        integer :: is

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi, weights)
        eta = exp(eta)

        ! calculate aggregate labor supply
        LL = sum(eta*weights)

        ! initialize grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize consumption policy function
        do is = 1, ns
            c(:, is) = 0.04d0*a(:) + eta(is)
            c_new(:, is) = c(:, is)
        enddo

        ! get an initial guess for the distribution function
        phi = 1d0/dble((NA+1)*NS)

    end subroutine
    
    

    subroutine value_function()

        use toolbox

        implicit none
        integer :: iter, ia, is, is_p
        real*8 :: V_new(0:NA, NS), aplus, con_lev
        integer :: ial(0:NA, NS), iar(0:NA, NS)
        real*8 :: varphi(0:NA, NS)

        ! get the interpolation shares and points
        do ia = 0, NA
            do is = 1, NS

                ! calculate where this guy would go
                aplus = max(((1d0+rn)*a(ia)+w*eta(is)-c(ia,is)), a_l)

                ! determine the gridpoints in between this decision lies and share
                call linint_Grow(aplus, a_l, a_u, a_grow, NA, &
                    ial(ia, is), iar(ia, is), varphi(ia, is))
            enddo
        enddo

        ! initialize value function
        V = c

        ! iterate until the value function converges
        do iter = 1, itermax

            ! calculate new value function
            do ia = 0, NA
                do is = 1, NS

                    ! interpolate over all future states
                    V_new(ia, is) = 0d0
                    do is_p = 1, NS
                        V_new(ia, is) = V_new(ia, is) + pi(is, is_p)* &
                            (varphi(ia,is)*V(ial(ia,is),is_p) + &
                            (1d0-varphi(ia,is))*V(iar(ia,is),is_p))**egam
                    enddo

                    ! add current utility
                    V_new(ia, is) = (c(ia,is)**egam + beta*V_new(ia, is))**(1d0/egam)
                enddo
            enddo

            con_lev = maxval(abs(V_new(:, :) - V(:, :))/max(abs(V(:, :)), 1d-10))

            ! update value function
            V = V_new

            ! check for convergence
            if(con_lev < sig_in)then
                V = V**egam/egam
                return
            endif
        enddo

    end subroutine


  
end program
