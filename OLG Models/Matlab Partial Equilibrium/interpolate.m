function [RHS EV]=interpolate(ij, cons, labor, V)
global a JJ JR NP NS NA egam nu beta pi gamma theta eta  eff pen w r psi
        
    RHS=zeros (JJ, NA+1, NP, NS); EV=zeros(JJ, NA+1, NP, NS);
    
       for ia=1:NA+1
            for ip = 1:NP
                for is=1:NS
                    
                       %RHS=zeros(ij, ia, ip, is);
                       %EV=zeros(ij, ia, ip, is);
                        
                    for is_p=1:NS
                       chelp=max(cons(ij,ia, ip,is_p),1d-10);
                       lhelp=max(labor(ij,ia, ip, is_p),1d-10);
                       %if pi(is, is_p)== 0
                       %       RHS(ij, ia, ip, is)=0;
                       %       fprintf('ij=%i, ia=%i, ip=%i, is=%i, is_p=%i,  RHS(ij, ia, ip, is)=%i, chelp=%c, lhelp=%i. \n', ij, ia, ip, is, is_p,RHS(ij, ia, ip, is),chelp, lhelp)
                       %else
                       RHS(ij, ia, ip, is)=RHS(ij, ia, ip, is)+ pi(is, is_p)*margu(chelp,lhelp);
                      % end
                       if  isnan(RHS(ij, ia, ip, is))
                            %disp(['ij= ', ij, ', ia', ia, ',ip=', ip, ',is=', is, ',is_p',is_p ])        
                            fprintf('ij=%i, ia=%i, ip=%i, is=%i, is_p=%i,  RHS(ij, ia, ip, is)=%i, chelp=%c, lhelp=%i. \n', ij, ia, ip, is, is_p,RHS(ij, ia, ip, is),chelp, lhelp)
                             error('gives me back Nan')
                            RHS(ij, ia, ip, is)=inf;
                       end
                       
                       
                       EV(ij, ia, ip, is)=EV(ij,ia, ip, is) +pi(is, is_p)*V(ij, ia, ip, is_p);
                       % margu=nu*(chelp^nu*(1-lhelp)^(1-nu))^(1-1/gamma)/chelp
                    end
                
                
                    RHS(ij, ia, ip, is)=(beta*psi(ij)*(1d0+r)*RHS(ij, ia, ip, is))^(-gamma);
                    %fprintf('psi(ij)=%i \n', psi(ij))
                    EV(ij, ia, ip, is)=(egam*EV(ij, ia, ip, is))^(1/egam);
                    
                      
                 end                  
            end         
       end
end



%GGG=ij=34; ia=33, ip=2, is=1, is_p=7,  RHS(ij, ia, ip, is)=NaN, chelp=1.000000e-10, lhelp=1