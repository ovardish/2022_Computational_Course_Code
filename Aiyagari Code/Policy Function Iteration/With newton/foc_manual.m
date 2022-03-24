function foc=foc_manual(x, COM)
       global          r s gridk mu
       i_com=COM(1);
       j_com=COM(2);
       r=COM(3);
       wage=COM(4);
       RHS_help=COM(5:end);
       
       %% block foc which puts two equations into fzero to solve;
          
       % calculate consumption
       a_plus = (1+r)*gridk(i_com) + s(j_com)*wage - x ;
       % calculate linear interpolation for future part of first order condition
       a_plus= min(max(a_plus, gridk(1)),gridk(end));
       % call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
       [ial, iar, varphi]=locate(a_plus, gridk);
       tomorrow = varphi*RHS_help(ial) +  (1d0-varphi)*RHS_help(iar);   
       
      % Uprime=margu(cons_com,lab_com);
       foc_int=x^(-mu) - tomorrow;
        
       %calculate first order condition for consumption       
       foc(1) = foc_int;
       % Get complementarity slackness condition !
end
