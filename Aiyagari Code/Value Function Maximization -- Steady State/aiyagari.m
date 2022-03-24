function [meank] = aiyagari(r)
% aiyagari2.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper

% r is interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta mu delta A alpha s N prob b probst kk gridk v tiny
%  
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
                                
                              
% form capital grid 
maxgridk = 25;                          % maximum value of capital grid  
mingridk = -phi;                        % borrowing constraint
incgridk = 0.09;                        % size of capital grid increments
gridk    = mingridk:incgridk:maxgridk;  % state of assets 
ngridk   = length(gridk);               % number of grid points

%  initialize some variables
test    = 10;
   
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum         
v   = zeros(N,ngridk);

%contemporaneous utility 
for j=1:N
    for i=1:ngridk  
        % constructing consumption for a productivity s(j), assets gridk
        % (i) and future assets gridk(l) where l=1, ngridk!
        cons = s(j)*wage + (1+r)*gridk(i) -gridk';  % cons dimension(Nx ngridk x ngridk)
        % Calculate the utility function      
        utilm(j,i,:) = (max(tiny,cons).^(1-mu)-1)./(1-mu);
    end
end
while test > 0.001;
    for j=1:N
        for i=1:ngridk              
            util=reshape(utilm(j,i,:),ngridk,1);            
            vint = util' + beta*prob(j,:)*v;         
           [tv(j,i),tdecis(j,i)] = max(vint);
        end
    end
    test=max(max(abs(tv-v)));
    v=tv;
 end;
  
%-----------------------------------------------------------------------
   %   form transition matrix
   %   trans is the transition matrix from state at t (row)
   %   to the state at t+1 (column) 
   %   The eigenvector associated with the unit eigenvalue
   %   of trans' is  the stationary distribution. 

   
% %% Method 1  start here
%        trans =spalloc(N*ngridk,N*ngridk,N*ngridk*2*N);
% 
%    for i = 1:ngridk      
%       for j = 1:N
%          trans((i-1)*N+j,(tdecis(j,i)-1)*N+1:tdecis(j,i)*N) = prob(j,:)';
%       end
%    end
%    
%    trans=trans';
%    
%    probst = (1/(N*ngridk))*ones(N*ngridk,1);
%    
%    test=1;
%    while test > 10^(-5);
%       probst1 = trans*probst;  
%       test = max(abs(probst1-probst));
%       probst = probst1;
%    end;
% %% Method 1  end here
  
   
  
%% Method 2  start here %%
lambda     = ones(N, ngridk)/(N*ngridk);
lambda_new = zeros(N,ngridk);
err=1;
   while err>tiny
   for j=1:N                % Current productivity state (in t)
       for i=1:ngridk       % Current asset state (in t)
           for j_p = 1:N    % Future productivity state (in t+1)       
           %   Asset state (in t+1)           
           lambda_new(j_p,tdecis(j,i))=  lambda_new(j_p,tdecis(j,i))+lambda(j, i)*prob(j,j_p); 
           end
       end
   end
   err=max(abs(lambda_new(:)-lambda(:)));
   lambda=lambda_new;
   lambda_new = zeros(N,ngridk);
   err;
 end
probst=lambda(:);
%% Method 2 end here %% 



%   vectorize the decision rule to be conformable with probst
%   calculate new aggregate capital stock  meanK
kk=gridk(tdecis(:));
meank=kk*probst;
   