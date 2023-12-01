function foc=foc_manual(x, COM)
       global a JJ JR NP NS NA egam nu beta pi gamma theta eta  eff pen w r psi
       ij_com=COM(1);
       ia_com=COM(2);
       ip_com=COM(3);
       is_com=COM(4);
       RHS_help=COM(5:end);
       
       %% block foc which puts two equations into fzero to solve;
       %x(1)=x_in;            
       %calculate the wage rate
       wage = w*eff(ij_com)*theta(ip_com)*eta(is_com);      
       % calculate available resources
       available = (1d0+r)*a(ia_com) + pen(ij_com);        
       % determine labor
       if (ij_com < JR)             
       lab_com = min( max( (1d0-nu)*(x(1) - available)/wage + nu, 1d-10) , 1d0);
       else
       lab_com = 0d0;
       end

       % calculate consumption
       cons_com = max( (available + wage*lab_com - x(1)) , 1d-10);
       % calculate linear interpolation for future part of first order condition
       a_plus= min(max(x, a(1)),a(end));
       % call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
       [ial, iar, varphi]=locate(a_plus, a);
       tomorrow = varphi*RHS_help(ial) +  (1d0-varphi)*RHS_help(iar);   
       
       Uprime=margu(cons_com,lab_com);
       foc_int=Uprime^(-gamma) - tomorrow;
        
       %calculate first order condition for consumption       
       foc(1) = foc_int;
       % Get complementarity slackness condition !
end
