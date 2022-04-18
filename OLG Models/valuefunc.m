function   [Value]= valuefunc(a_plus, cons, lab, ij, ip, is, EV)
       global a JJ JR NP NS NA egam nu beta pi gamma theta eta  eff pen w r psi
        %integer, intent(in) :: ij, ip, is
        %real*8, intent(in) :: a_plus, cons, lab
        %real*8 :: valuefunc, varphi, c_help, l_help
        %integer :: ial, iar
        % check whether consumption or leisure are too small
        c_help = max(cons, 1d-10);
        l_help = min(max(lab, 0d0),1d0-1d-10)  ; 
%%
        % get tomorrows utility
        [ial, iar, varphi]=locate(a_plus, a); % Linear interpolant !  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        
        %! calculate tomorrow's part of the value function
        Value = 0d0;
        if  ij < JJ
            Value = max(varphi*EV(ij+1, ial, ip, is) +...
                (1d0-varphi)*EV(ij+1, iar, ip, is), 1d-10)^egam/egam ;  % Varphi is a interpolation weights ! I need to put here splines ! Probably, but for now we live it herre !!!
        end
        
        % add todays part and discount
        Value = (c_help.^nu*(1d0-l_help).^(1d0-nu)).^egam/egam + beta*psi(ij+1)*Value;
    end 
