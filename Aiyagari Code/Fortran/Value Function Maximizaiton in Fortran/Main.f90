program Main
    
    use toolbox
    use globals
    !use ModuleSaving
    
    implicit none
    integer :: solution_method , iter
    real*8 ::r_new    , minrate, maxrate, KK_0
    
    
    minrate     = -delta;
    maxrate     = (1-beta)/(beta); 
    

    
     do  iter=1, itermax
         
     ! set initial guess of the interest rate
      r = 0.5d0*(minrate+maxrate)   

     ! initialize variables
     if (iter==1)      call initialize() 
    
     KK_0      = ((r+delta)/(alpha*LL**(1-alpha)))**(1/(alpha-1));     ! implied guess for capital stock   
     w=(1d0-alpha)*((r+delta)/(alpha))**(alpha/(alpha-1d0))       
    
     ! Here solve the household problem   
     ! Value function maximizaiton with fminsearch      
     call minimize
     
     ! calculate the invariant distribution directly
      call get_distribution()
     
     ! get aggregate assets
      AA = sum(a*sum(phi, 2)   )   
      KK=AA

   
 
     
    
    
    print*, '**  r', r
    print*, '**  K_demand', KK_0
    print*, '**  K_supply', KK
    ! Main Aiyagari loop nested in the capital bisection loop
    !update guesses
    ! find optimal policies
    if ( abs(KK-KK_0) <0.00001) then        
             call Save_policies
             print*, 'Found Steady State'
        exit
    endif
      
    !! updating the interest rate limits
    if ( KK >  KK_0   ) then    ! If the supply of capital stock is higher than demand, lower the interest rate  
        
          maxrate=r;
    else
          minrate=r;             ! If the supply is lower than demand, than increase the interest rate
    endif
    
    enddo
        
  
          
   contains
         
    ! For initializing variables
    subroutine initialize()
        use toolbox
        implicit none
        integer :: is        
        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta_p, pi_p, weights_p)
        eta_p = exp(eta_p)
        if (NS_p==NS) then
        eta(1:NS_p)=eta_p(1:NS_p)
        pi(1:NS_p, 1:NS_p) = pi_p(1:NS_p, 1:NS_p)
        weights(1:NS_p)=weights_p(1:NS_p)
     else
        eta(NS-NS_p+1:NS) = eta_p
        eta(1) =0.0d0
       pi(NS-NS_p+1:NS, NS-NS_p+1:NS) = pi_p
       pi(1,NS-NS_p+1:NS)=weights_p
       pi(:,1)=0.1d0
         
        do is=1,NS
            pi(is,:)=pi(is,:)/sum(pi(is,:))
            print*, sum(pi(is, :))
        enddo   
     endif
        ! calculate aggregate labor supply
        LL = sum(eta_p*weights_p/sum(pi(1,:)))
        ! initialize grid with equidistant points
         ! if ( grid_method==1)  then
        !         call grid_Cons_Equi(a, a_l, a_u)          ! Equidistant Grid
        ! else
                 call grid_Cons_Grow(a, a_l, a_u, a_grow)
       !  endif
        
         
        ! initialize policy function
        do is = 1, ns            
            c_pol(:, is) = r*a(:) + eta(is)*w
            V(:, is) = c_pol(:,is)**(egam)/(egam)
        enddo
             EV=(egam*V)**(1.0d0/egam)

        ! get an initial guess for the distribution function
        phi = 1d0/dble((NA+1)*NS)
            
    end subroutine

          
       
    
   subroutine Minimize


    implicit none
    integer :: ia, iter , is  , iter_howard, howard_ind
    
    howard_ind=1
    ! start timer
    call tic()

    ! interpolate coefficients
        

 
    ! iterate until value function converges
    do iter = 1, itermax

        ! set a = 0 manually
        c_pol(0,1) = 0d0
        V_new(0,1) = -1d10

        ! calculate optimal decision for every gridpoint
    !$OMP PARALLEL DO DEFAULT(SHARED)
        do ia = 0, NA
           do is=1, NS               
            ! initialize starting value and communicate resources  given the initial guess of c_pol
            x_in = a(ia)*(1+r) +w*eta(is) - c_pol(ia,is)
            ! to pass the function utility and fminsearch, you need to set up communication variables
            ia_com = ia
            is_com = is
            
            ! find the maximm of the function UTILITY given the state space (ia, is); Here I am just using a smart guess, that the chosen capital cannot be more than the available resources in hand, assuming c_pol=0
            call fminsearch(x_in, fret, a_l, a(ia)*(1+r) +w*eta(is), utility)

           !  stop
            ! get optimal consumption and value function
            if (x_in<0d0)      print*, x_in
            c_pol(ia, is) = a(ia)*(1+r) +w*eta(is) - x_in
            V_new(ia, is) = -fret
            a_pol(ia, is) = x_in
            
            
           enddo
        enddo

        
     ! This is to parallelize, we will use it in the future
     !$OMP END PARALLEL DO
     ! interpolate coefficients
       !$OMP PARALLEL DO DEFAULT(SHARED)
        EV=0.0d0
        do ia=0, NA            !        
         do is=1,NS 
         EV(ia, is) =EV(ia, is)+ sum(pi(is, :) * V_new(ia, :)) 
         enddo
       enddo
      !$OMP END PARALLEL DO
      EV=(egam*EV)**(1.0d0/egam)      

      
        
     ! get convergence level
      con_lev = maxval(abs(V_new(:,:) - V(:,:))/max(abs(V(:,:)), 1d-10))
      write(*,'(i5,2x,f20.7)')iter, con_lev
      
     ! check for convergence
     if(con_lev < sig_in*100)then
     print*, ' c o n v e r g e d  '         
     exit   
     endif

        
   ! Hoover algorithm  
   if (howard_ind==1)     then
   do iter_howard = 1, 350 ! Howard
       ! Now we will iterate over the policy function
      !$OMP PARALLEL DO DEFAULT(SHARED)
        do ia = 0,NA
           do is = 1, NS
               ia_com=ia
               is_com=is
                V_new(ia,is) = - utility(a_pol(ia,is))
            enddo
        enddo
       !$OMP END PARALLEL DO
      ! interpolate coefficients
          EV=0.0d0
       !$OMP PARALLEL DO DEFAULT(SHARED)
       do ia=0, NA            !        
         do is=1,NS 
         EV(ia, is) =EV(ia, is)+ sum(pi(is, :) * V_new(ia, :)) ! max(spline_eval(x_in, coeff_v(:,is_p), a_l, a_u), 1d-10)**egam/egam
       
         enddo
      enddo
      !$OMP END PARALLEL DO
      EV=(egam*EV)**(1.0d0/egam)
       enddo
   endif  ! hoover index   
   
   V=V_new

   enddo
    
   end subroutine    


    
    
      subroutine save_policies
       open(1,  file = 'Output/V_new.txt',      status = 'unknown')
        open(2,  file = 'Output/c_pol.txt',      status = 'unknown')
        open(3,  file = 'Output/a.txt',      status = 'unknown')
        open(4,  file = 'Output/a_pol.txt',      status = 'unknown')
        open(5,  file = 'Output/phi.txt',      status = 'unknown')
        open(6,  file = 'Output/EV.txt',      status = 'unknown')
         
        write(1, *) V_new
        write(2, *) c_pol
        write(3, *) a
        write(4, *) a_pol
        write(5, *) phi
        write(6, *) EV
          
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
      
      end subroutine
    
   
    
end program
