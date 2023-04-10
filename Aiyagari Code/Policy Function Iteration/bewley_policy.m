%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aiyagari's model with Policy Function Iterations
% contact for typos: Oliko.Vardishvili@uci.edu
% Codes based on S&L2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
global beta mu A delta alpha s b N prob gridk probst kpol c1

%
%  set parameter values
%
mu     = 3;               % risk aversion              
beta   = 0.96;            % subjective discount factor 
delta  = 0.08;            % depreciation
A      = 1;            % production technology
alpha  = 0.36;            % capital's share of income

% approximate labor endowment shocks with seven states Markov chain
% log(s_t) = rho*log(s_t-1)+sigmaint*sqrt(1-rho^2)*error_t

N        = 7;             % number of discretized states
rho      = 0.2;           % first-order autoregressive coefficient
sigmaint = 0.4;           % intermediate value to calculate sigma
sigma    = sigmaint*sqrt(1-rho^2); % standard deviation of error_t

% prob is transition matrix of the Markov chain
% logs is the discretized states of log labor earnings
% invdist is the invariant distribution of Markov chain
% alambda and asigma are respectively the theoretical 
% rho and standard deviation of log labor income in Markov chain

[prob,logs,invdist,alambda,asigma]=markovappr(rho,sigma,3,N);

s = exp(logs);
labor = s*invdist;

b =2;           %BORROWING CONSTRAINT

%setting limits for interest rate search
minrate = -0.01;
maxrate = (1-beta)/beta;

%setting test values
testr=1;
testk=1;
% we find the optimal interest rate by bisection
tic
while (testk>0.0001)
    r0 = 0.5*(maxrate-minrate)+minrate;
    k0=((r0+delta)/(alpha*A*labor^(1-alpha)))^(1/(alpha-1));    
    k1= aiyagari_policy_newton(r0); %calculating asset demand    
    r1=alpha*A*max(0.001,k1)^(alpha-1)*labor^(1-alpha)-delta;
    [k0 k1 r0 r1];
    testk=abs(r1-r0);    
    %updating the interest rate
    if k1>k0
        maxrate=r0;
    else
        minrate=r0;
    end
    
    disp(['k0 = ',num2str(k0),', k1 = ',num2str(k1),', r0 = ',num2str(r0),', r1 = ',num2str(r1)])

end
steady.time = toc 
fprintf('Steady State solved, it took me %3.4f seconds. \n\n',[steady.time])

ngridk=length(gridk);
dist=reshape(probst,N,ngridk);
wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));
dist=reshape(probst,N,ngridk);


%% Formating Figures
% Fontsize
FS=12;
% Linewidth
WD=2.5;
figure(1)
plot(gridk,[gridk' kpol(:,1) kpol(:,4) kpol(:,7)], '-.','linewidth', WD)
ylabel('next period asset holdings')
xlabel('current asset holdings')
set(gca,'XGrid','off','YGrid','on','Fontsize',FS)
set(gca,'TickLabelInterpreter','LaTex')
xlim([-3 gridk(ngridk)])

figure(2)
plot(gridk,[dist(4,:); dist(5,:); sum(dist); ],  '-','linewidth', WD)
ylabel('distribution')
xlabel('asset holdings')
set(gca,'XGrid','off','YGrid','on','Fontsize',FS)
set(gca,'TickLabelInterpreter','LaTex')
xlim([-3 gridk(ngridk)])

figure(3)
plot(gridk,[c1(:,1) c1(:,4) c1(:,7)], '-.','linewidth', WD)
ylabel('consumption')
xlabel('asset holdings')
set(gca,'XGrid','off','YGrid','on','Fontsize',FS)
set(gca,'TickLabelInterpreter','LaTex')
xlim([-3 gridk(ngridk)])



