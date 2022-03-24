
function [meank] = aiyagari_policy_fehr(r)
% aiyagari2.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper
% r is interest rate


global beta mu delta A alpha s N prob b probst kk gridk kpol c1
%   write wage as a function of interest rate using Cobb-Douglas
%   production function

wage = (1-alpha)*(A*(alpha/(r+delta))^alpha)^(1/(1-alpha));

if r<=0
   phi = b;
else
   phi  = min(b, wage*s(1)/r);           
end

% -phi is borrowing limit, b is adhoc
% the second term is natural limit

%--------------------------------------------------------------------                               
                                
                              
%   form capital grid

%   
   
maxgridk = 25;                     % maximum value of capital grid  
mingridk = -phi;                   % borrowing constraint
incgridk = 0.09;                    % size of capital grid increments
gridk    = mingridk:incgridk:maxgridk;   % state of assets 
ngridk   = length(gridk);            % number of grid points

%  initialize some variables
%
c0       = zeros(ngridk,N);  %consumption function   

test    = 10;
   
%  iterate on the Euler equations and get the decision 

for j=1:N
    for i=1:ngridk              
        c0(i,j) = max(0.01,s(j)*wage + r*gridk(i)); % initial guss for consumption
        wealth(i,j)=s(j)*wage + (1+r)*gridk(i);     % defining wealth, i.e., for each state, what is available resources? 
    end
end
c1=c0;

while test > 0.001
       
    RHS=beta*(1+r).*c0.^(-mu)*prob';

    for j=1:N
       for i=1:ngridk
             
              x=c0(i,j);
              COM=[i, j, r,wage, RHS(  :,j)'];
              options = optimoptions('fsolve','Display','off');
              f=@(x) foc_manual(x, COM); 
              [X, fval]=fsolve(f, x,options);
                                  
              if  X > (gridk(i)*(1+r) +s(j)*wage - gridk(1))     ;          
                  X =  (gridk(i)*(1+r) +s(j)*wage - gridk(1));
              end
              
              c1(i,j)=X;
              a1(i,j)=gridk(i)*(1+r) +s(j)*wage -   c1(i,j);
                   
        end
    end
    
    test=max(max(abs(c1-c0)))
    c0=0.1.*c0+0.9.*c1;
    
end
 
 
 
 

%-----------------------------------------------------------------------
   %   form transition matrix
   %  now the states are continuos, see hand-out for interpretation
   trans(N*ngridk,N*ngridk)=0;
   
  for j = 1:N
       for i = 1:ngridk      
            [xx indopt(j,i)]=min(abs(kpol(i,j)-gridk));
            if kpol(i,j)>gridk(indopt(j,i)) | kpol(i,j)==gridk(1)
                indopt(j,i)=indopt(j,i)+1;
            end
       end
  end
  
  for j = 1:N
      for i = 1:ngridk      
          if i==1
               if (kpol(i,j)==gridk(1))
                   trans((i-1)*N+j,1:N) = prob(j,:)';
               else                         
                   trans((i-1)*N+j,(indopt(j,i)-1)*N+1:indopt(j,i)*N) = prob(j,:)';
               end             
          else
             if (kpol(i-1,j)==gridk(1) & kpol(i,j)==gridk(1))
                 trans((i-1)*N+j,1:N) = prob(j,:)';
             else
                 if indopt(j,i)==indopt(j,i-1)
                    trans((i-1)*N+j,(indopt(j,i)-1)*N+1:indopt(j,i)*N) = prob(j,:)';
                 else
                     for l=indopt(j,i-1):indopt(j,i)
                        trans((i-1)*N+j,(l-1)*N+1:l*N) = prob(j,:)'*(min(gridk(l),kpol(i,j))-max(kpol(i-1,j),gridk(l-1)))/(kpol(i,j)-kpol(i-1,j));
                     end
                 end
             end
        end
      end
  end
  
   trans=trans';
   probst = (1/(N*ngridk))*ones(N*ngridk,1);
   
   test=1;
   while test > 10^(-7);
      probst1 = trans*probst;   
      test = max(abs(probst1-probst));
      probst = probst1;
   end;
   %
   %   vectorize the decision rule to be conformable with probst
   %   calculate new aggregate capital stock  meanK
   %   

   kpols=kpol;
   %using average capital accumulation for the intervals
   kpols(2:ngridk,:)=0.5*kpol(1:ngridk-1,:)+0.5*kpol(2:ngridk,:);
   kpols=kpols';
   kk=kpols(:);
   meank=kk'*probst;
   
