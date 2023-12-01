clear all
close all
clc

global a JJ JR NP NS NA egam nu beta pi gamma theta eta  eff pen w r psi

%% Define Numerical variables;

% number of years the household retires
JR=45; 
% number of years the household lives
JJ=80;
% number of persistent shock process values
NP=2; %Fixed effect theta
% number of transitory shock process values
NS=7;
% number of points on the asset grid
NA=200;    

% household preference parameters HH utility: [c_j^v(1-l_j)^(1-v)]^(1-1/gamma)/(1-1/gamma)
gamma  = 0.50d0;
egam   = 1d0 - 1d0/gamma;
nu     = 0.335d0;
beta   = 0.98d0;

% household risk process
sigma_theta = 0.242d0;  % Fixed effect
sigma_eps   = 0.022d0;  % Idiosyncratic shock: eta_j=rho eta_j-1+epsilon_j
rho         = 0.985d0;

% Size of the asset grid, note we do not do equidistant grid: a_v=a_l+(a_u-a_l) * ((1+a_grow)^v-1)/((1+a_grow)^n-1)  for v=0,1, ..., n
% Note: (a_{v+2}-a_{v+1})/(a_{v}-a_{v-1})-1=a_grow
a_l    = 0.0d0;   
a_u    = 200d0;
a_grow = 0.05d0;
    
% net prices
%r, w

                       
% transfer payments (old-age) and survival probabilities
% pen(JJ), psi(JJ+1)

%cohort aggregate variables
c_coh=zeros(JJ,1); n_coh=zeros(JJ,1); l_coh=zeros(JJ,1); h_coh=zeros(JJ,1); a_coh=zeros(JJ,1); v_coh=zeros(JJ,1);


% the permanent shock process
dist_theta=ones(NP,1); theta=zeros(NP,1);

% the transitory shock process
pi=zeros(NS, NS);
eta=zeros(NS,1);   % Markov Chain probability matrix & states
is_initial = 4;    %     
   
   
% demographic and other model parameters
eff=zeros(JJ,1)    
    
    
% individual variables
a=zeros(NA+1); aplus=ones(JJ, NA+1, NP, NS);   % a asset today, aplus _ calculated future assets!
c=zeros(JJ, NA+1, NP, NS); l=zeros(JJ, NA+1, NP, NS);
phi=zeros(JJ, NA+1, NP, NS); V=zeros(JJ, NA+1, NP, NS);

% numerical variables
RHS=zeros (JJ, NA+1, NP, NS); 
EV=zeros(JJ, NA+1, NP, NS);
%ij_com, ia_com, ip_com, is_com
%cons_com, lab_com
 
  
  
%% subroutine initialize        
r = 0.04d0;
w = 1.0d0;

% set survival probabilities
psi = [1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, ...
                0.99906d0, 0.99908d0, 0.99906d0, 0.99907d0, 0.99901d0, ...
                0.99899d0, 0.99896d0, 0.99893d0, 0.99890d0, 0.99887d0, ...
                0.99886d0, 0.99878d0, 0.99871d0, 0.99862d0, 0.99853d0, ...
                0.99841d0, 0.99835d0, 0.99819d0, 0.99801d0, 0.99785d0, ...
                0.99757d0, 0.99735d0, 0.99701d0, 0.99676d0, 0.99650d0, ...
                0.99614d0, 0.99581d0, 0.99555d0, 0.99503d0, 0.99471d0, ...
                0.99435d0, 0.99393d0, 0.99343d0, 0.99294d0, 0.99237d0, ...
                0.99190d0, 0.99137d0, 0.99085d0, 0.99000d0, 0.98871d0, ...
                0.98871d0, 0.98721d0, 0.98612d0, 0.98462d0, 0.98376d0, ...
                0.98226d0, 0.98062d0, 0.97908d0, 0.97682d0, 0.97514d0, ...
                0.97250d0, 0.96925d0, 0.96710d0, 0.96330d0, 0.95965d0, ...
                0.95619d0, 0.95115d0, 0.94677d0, 0.93987d0, 0.93445d0, ...
                0.92717d0, 0.91872d0, 0.91006d0, 0.90036d0, 0.88744d0, ...
                0.87539d0, 0.85936d0, 0.84996d0, 0.82889d0, 0.81469d0, ...
                0.79705d0, 0.78081d0, 0.76174d0, 0.74195d0, 0.72155d0, ...
                0.00000];
            
 psi=psi';

% initialize age earnings process
eff(1:JR-1) =[1.0000d0, 1.0719d0, 1.1438d0, 1.2158d0, 1.2842d0, 1.3527d0, ...
              1.4212d0, 1.4897d0, 1.5582d0, 1.6267d0, 1.6952d0, 1.7217d0, ...
              1.7438d0, 1.7748d0, 1.8014d0, 1.8279d0, 1.8545d0, 1.8810d0, ...
              1.9075d0, 1.9341d0, 1.9606d0, 1.9623d0, 1.9640d0, 1.9658d0, ...
              1.9675d0, 1.9692d0, 1.9709d0, 1.9726d0, 1.9743d0, 1.9760d0, ...
              1.9777d0, 1.9700d0, 1.9623d0, 1.9546d0, 1.9469d0, 1.9392d0, ...
              1.9315d0, 1.9238d0, 1.9161d0, 1.9084d0, 1.9007d0, 1.8354d0, ...
              1.7701d0, 1.7048d0];
eff(JR:JJ) = 0d0;

% old-age transfers
pen = 0d0;
pen(JR:JJ) = 0.5d0*sum(eff)/(JR-1)*0.33d0;

% initialize fixed effect
dist_theta = 1d0/NP.*dist_theta;
theta(1)   = exp(-sqrt(sigma_theta));  % one standard deviation
theta(2)   = exp( sqrt(sigma_theta));

% Idiosyncratick shocks 
%[eta,pi] = tauchen(NS,0,rho,sigma_eps,4);
%eta=exp(eta);
eta =  [0.121781321698703;       0.245689797353754;       0.495671057611551;
   1.00000000000000;        2.01746699679948;        4.07017308317513;
   8.21143986656742;]      

pi=  [0.955835359818732 7.222937227849423E-003 5.458138963110441E-005 4.124538259277446E-007 3.116779541015728E-009 2.355249023437597E-011 1.779785156250089E-013; ...
      4.333762336709654E-002  0.95610826676688 1.444752427100256E-002  1.637535192319363E-004  1.649862408691447E-006 1.558407568359426E-008  1.413149414062559E-010; ... 
      8.187208444665662E-004 3.611881067750639E-002  0.956272029636458       2.167252383863909E-002 3.275132722009331E-004  4.124656021728618E-006  4.675169311523592E-008; ...
      8.249076518554891E-006  5.458450641064543E-004  2.889669845151879E-002 0.956326620376606      2.889669845151879E-002  5.458450641064543E-004 8.249076518554891E-006; ...
      4.675169311523592E-008  4.124656021728618E-006 3.275132722009331E-004  2.167252383863909E-002  0.956272029636458 3.611881067750639E-002  8.187208444665662E-004;...
      1.413149414062559E-010  1.558407568359427E-008  1.649862408691447E-006  1.637535192319363E-004 1.444752427100256E-002  0.956108266766888       4.333762336709654E-002; ...
      1.779785156250089E-013  2.355249023437597E-011  3.116779541015728E-009  4.124538259277446E-007  5.458138963110441E-005  7.222937227849423E-003  0.955835359818732;]

pi=pi';
 
  
  


%% Grid for assets 
% Discretize Grid
%call grid_Cons_Grow(a, a_l, a_u, a_grow)
%n = size(a, 1)-1
% check for left <= right
%if(a_l <= a_u)
%error('grid_Cons_Grow: a_l must be less then a_u')
%end
%if(growth <= 0d0)
%error('grid_Cons_Grow: growth rate must be greater than zero')
%end

% calculate factor
%h = (a_r-a_l)/((1+a_grow)**n-1d0)
     
% calculate grid value
%for kk=1:n+1
%a = h*((1d0+growth)**j=0,n)]-1d0)+left

% Discretized asset grid
a_lb=log(1);
a_ub=log(a_u);
a=linspace(a_lb, a_ub, NA+1);
a=exp(a)-1;
x=zeros(2,1); % first entry for assets, seccond entry for slack variable.

%% SOLVE HOUSEHOLD! MOST IMPORTANT PART
% get decision in j=JJ;
  for ia=1:NA+1
       aplus(JJ, ia, :, :) = 0d0;
       c(JJ, ia, :, :) = (1d0+r)*a(ia) + pen(JJ);
       l(JJ, ia, :, :) = 0d0;
       V(JJ, ia, :, :) = valuefunc(0d0, c(JJ, ia, 1, 1), l(JJ, ia, 1, 1), JJ, 1, 1, EV);
  end
[RHS, EV]=interpolate(JJ, c, l, V);


EV(JJ,:,1,1,1,1)
%% Should i do it here, or even down?   
   for ij=JJ-1:-1:1
           if ij<JR
               ip_max=NP;
               is_max=NS;
           else
               ip_max=1;
               is_max=1;
           end
           


         
       for ia=1:NA+1
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%% FIRST I DO NOT KNOW IF IT IS RIGHT
           if ij>=JR && ia ==0 && pen(ij) <=1d-10 
               aplus(ij, ia, :, :)=0d0;
               c(ij,ia, :,:)=0d0;
               l(ij,ia, :,:)=0d0;
               V(ij,ia, :,:)=valuefunc(0d0,0d0, 0d0, ij, 1, 1);
               %cycle
           end
           
   
           
           for ip=1:ip_max
               for is=1:is_max
               x_in=aplus(ij+1, ia, ip, is);
               COM=[ij, ia, ip, is, RHS(ij+1, :, ip, is)];
               options = optimoptions('fsolve','Display','off');
               %f=@(x) focs(x, COM); 
               %[X, fval]=fsolve(f, [x_in 0.5],options);
               f=@(x) foc_manual(x, COM); 
               %[X, fval]=fzero(f, x_in);
               [X, fval]=fsolve(f, x_in,options);
               
               %if isnan(X)
               %    error
               %end
               
               if    X<0;
                     X=0;
                   %   wage = w*eff(ij)*theta(ip)*eta(is);      
                      %calculate available resources
                    %  available = (1d0+r)*a(ia) + pen(ij);        
                     % determine labor
                    %  if ij < JR             
                    %  l = min( max( nu-(1d0-nu)*(X - available)/wage, 1d-10) , 1d0);
                    %  else
                    %  l= 0d0;
                     % end
                     % calculate consumption
                    % c(ij,ia, ip, is) = max( (available + wage*l - X) , 1d-10);
                end
               
               
               %% so if we do not have complementarity slackness condition, we can easily check the bc violation x(1)<0, and manually rewrite it;  
               
               
               
               %% Store policy functions
               wage = w*eff(ij)*theta(ip)*eta(is); 
               aplus(ij,ia, ip, is)=X;
               
               available = (1d0+r)*a(ia) + pen(ij);        
               if ij<JR
               l(ij,ia, ip, is)= min( max( (1d0-nu)*(X - available)/wage + nu, 1d-10) , 1d0);
               else
               l(ij,ia, ip, is)=0;
               end
               c(ij,ia, ip, is)=max( (available + wage*l(ij,ia, ip, is) - X) , 1d-10);           
               V(ij,ia, ip, is)=valuefunc(X(1),c(ij,ia, ip, is), l(ij,ia, ip, is), ij, ip, is, EV);
             end
         end
           
         
         
         % Store decisions in retirement age
         if ij>= JR
             aplus(ij, ia, :,:)=aplus(ij,ia,1,1);
             c(ij,ia,:,:)=c(ij,ia,1,1);
             l(ij,ia,:,:)=l(ij, ia, 1, 1);
             V(ij, ia, :,:)=V(ij, ia, 1,1);
         end
           
       end
       
     
       [RHS, EV]=interpolate(ij, c, l, V);
      if ij==19
          plot(a,c(ij,:,ip, is), 'blue',a,c(ij+1,:,ip, is), 'red',a,c(ij+2,:,ip, is), 'red', a,c(ij+3,:,ip, is), 'yellow' )
          pause
       end
       fprintf('Age %i %i is done!\n', ij, is);
   
   end 
  
   
   
   
   
   
   
%% Get Distribution
   is_initial=(NS+1)/2;
   phi(:,:,:,:)=0d0;
   for ip=1:NP
       phi(1,1,ip, is_initial)=dist_theta(ip);
   end
   
% succesively compute distribution over ages
for ij=2:JJ
    for ia=1:NA+1
        for ip=1:NP
            for is=1:NS
                 [ ial, iar,  varphi] = locate(aplus(ij-1,ia,ip,is), a);
                % ial=min(ial,NA);
                % iar=min(iar, NA);
                 %varphi=min(varphi,1d0);
                 for is_p=1:NS
                     phi(ij,ial,ip,is_p)=phi(ij,ial,ip,is_p)+pi(is,is_p)*varphi*phi(ij-1, ia, ip, is);
                     
                     phi(ij,iar,ip,is_p)=phi(ij,iar,ip,is_p)+pi(is,is_p)*(1d0-varphi)*phi(ij-1, ia, ip, is);

                 end
            end
        end
    end
end

%% Get aggregation across cohorts;
for ij=1:JJ
    for ia=1:NA+1
       for ip=1:NP 
         for is=1:NS    
             c_coh(ij)=c_coh(ij)+c(ij,ia,ip,is)*phi(ij, ia, ip, is);
             n_coh(ij)=n_coh(ij)+l(ij, ia, ip, is)*phi(ij, ia, ip, is);
             l_coh(ij)=l_coh(ij)+eff(ij)*theta(ip)*eta(is)*l(ij, ia, ip, is)*phi(ij, ia, ip, is);
             h_coh(ij)=h_coh(ij)+eff(ij)*theta(ip)*eta(is)*phi(ij, ia, ip, is);
             a_coh(ij)=a_coh(ij)+a(ia)*phi(ij, ia, ip, is);
             v_coh(ij)=v_coh(ij)+V(ij, ia, ip, is)*phi(ij, ia, ip, is);             
         end
       end
    end 
end



%% Calculate cohort specific coefficients of variation
cv_c=zeros(JJ,1);
cv_n=zeros(JJ,1);
cv_l=zeros(JJ,1);
cv_h=zeros(JJ,1);
corr_hn=zeros(JJ,1); %! Coefficient of Variation

for ij=1:JJ
    for ia=1:NA+1
       for ip=1:NP 
         for is=1:NS    
             cv_c(ij)=cv_c(ij)+c(ij,ia,ip,is)^2*phi(ij, ia, ip, is);
             cv_n(ij)=cv_n(ij)+l(ij, ia, ip, is)^2*phi(ij, ia, ip, is);
             cv_l(ij)= cv_l(ij)+(eff(ij)*theta(ip)*eta(is)*l(ij, ia, ip, is))^2*phi(ij, ia, ip, is);
             cv_h(ij)=cv_h(ij)+(eff(ij)*theta(ip)*eta(is))^2*phi(ij, ia, ip, is);
             corr_hn(ij)=corr_hn(ij)+  eff(ij)*theta(ip)*eta(is)*l(ij,ia,ip,is)*phi(ij,ia,ip,is);           
         end
       end
    end 
end
corr_hn=(corr_hn-n_coh.*h_coh)./max(sqrt(cv_n-n_coh.^2).*sqrt(cv_h-h_coh.^2),1d-10);
cv_c=sqrt(cv_c-c_coh.^2)./c_coh;
cv_h=sqrt(cv_h-h_coh.^2)./max(h_coh, 1d-10);
cv_l=sqrt(cv_l-l_coh.^2)./max(l_coh,1d-10);
cv_n=sqrt(cv_n-n_coh.^2)./max(n_coh,1d-10);






% oliko=reshape(V(JJ, 200, :, 1), 2,[])
figure(1)
plot(1:JJ, c_coh, 'black', 1:JJ, n_coh, 'red', 1:JJ, l_coh+pen', 'blue', 'Linewidth', 3.25)
xlabel('$Age j$','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
%xlim([a_vec(1)-0.5 a_vec(Na)+0.5])
%ylim([a_vec(1)-0.5 max(a_pol)+0.5])

 
figure(1)
plot(1:JJ,c_coh,  1:JJ, n_coh,  1:JJ, l_coh+pen',  'Linewidth', 3.25)
title('c_coh, l_coh, h_coh','Interpreter','LaTeX','Fontsize',23)
xlabel('Age j','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
leg = legend('Consumption', 'Earnings', 'Hours');
set(leg,'Interpreter','LaTex','Fontsize',14,'Location','SouthEast')
legend boxoff

figure(2)
plot(1:JJ,v_coh, 'Linewidth', 3.25)
title('Value function over lifecycle','Interpreter','LaTeX','Fontsize',23)
xlabel('Age j','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
leg = legend('assets');
set(leg,'Interpreter','LaTex','Fontsize',14,'Location','SouthEast')
legend boxoff

figure(3)
plot(1:JJ,a_coh, 'Linewidth', 3.25)
title('Assets over lifecycle','Interpreter','LaTeX','Fontsize',23)
xlabel('Age j','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
leg = legend('assets');
set(leg,'Interpreter','LaTex','Fontsize',14,'Location','SouthEast')
legend boxoff

figure(4)
plot(1:JJ,cv_c,1:JJ,cv_l,1:JJ,cv_n, 'Linewidth', 3.25)
title('Coefficient of variations','Interpreter','LaTeX','Fontsize',23)
xlabel('Age j','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
leg = legend('cv_c', 'cv_l', 'cv_n');
set(leg,'Interpreter','LaTex','Fontsize',14,'Location','SouthEast')
legend boxoff
 
figure(5)
plot(1:JJ,cv_c,1:JJ,cv_l,1:JJ,cv_n,1:JJ,corr_hn, 'Linewidth', 3.25)
title('Variance Decomposition of Earnings','Interpreter','LaTeX','Fontsize',23)
xlabel('Age j','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
leg = legend('cv_c', 'cv_l', 'cv_n', 'corr_hn');
set(leg,'Interpreter','LaTex','Fontsize',14,'Location','SouthEast')
legend boxoff
 
 

 %% So now lets check how many HHs are actually facing 

 
 frac_bor=zeros(JJ,1);
 for ij =1:JJ-1
     for ia=1:NA+1
         for ip=1:NP
             for is=1:NS
                 if  aplus(ij, ia, ip, is)<1d-6
                     frac_bor(ij)=frac_bor(ij)+phi(ij,ia,ip,is);
                 end
             end
         end
     end
  end
 
 
frac_bor(JJ)=1;

figure(6)
plot(1:JJ,frac_bor,  'Linewidth', 3.25)
title('fraction of constrained HHs','Interpreter','LaTeX','Fontsize',23)
xlabel('$a$','Interpreter','LaTex','Fontsize',21)
set(gca,'XGrid','off','YGrid','on','Fontsize',17)
set(leg,'Interpreter','LaTex','Fontsize',14,'Location','SouthEast')
legend boxoff













%% Subroutine that checks for the maximum gridpoint used
iamax=zeros(JJ,1);
 for ij =1:JJ
     for ia=1:NA+1
         for ip=1:NP
             for is=1:NS
                 if  phi(ij, ia, ip, is)>1d-8
                     iamax(ij)=ia;
                 end
             end
         end
     end
  end