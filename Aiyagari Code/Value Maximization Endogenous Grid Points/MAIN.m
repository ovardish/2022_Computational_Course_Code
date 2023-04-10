%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving steady state of the Aiyagari model                                                  %%%
%% Using Endogeneous Grid Points                                                               %%%
%% Advanced Macro; Econ 526; Instructor Oliko Vardishvili                                      %%%
%% Contact for typos: oliko.vardishvili@uci.edu                                               %%%
%% Initial code was provided Lukas Nord, a Phd Candidate, European University Institute        %%%      
%% We are grateful to him : )                                                                  %%%                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The codes are meant to illustrate solution methods as easily
% accessible as possible and are not necessarily optimized for maximum
% speed. Suggestions to further improve the codes are very welcome.

clear;
close all;
clc;


%% Initialization

% algorithm
    param.tol_r           = 1e-6;          % tolerance for convergence of interest rate
    param.tol_pol         = 1e-12;         % tolerance for convergence of policy functions
  
% parameters
    param.gamma           = 2;               % risk aversion 
    param.beta            = 0.96;            % subjective discount factor 
    param.delta           = 0.08;            % depreciation
    param.A               = 1;               % aggregate productivity
    param.alpha           = 0.36;            % capital's share of income
    param.nkap            = 300;             % number of asset grid points
    param.b               = -2;              % exogenous borrowing limit
    param.B               = 30;             % upper bounnd on assets
    grid.k                = linspace(param.b,param.B,param.nkap); % equally spaced grid
    %grid.k                = param.b + linspace(0,1,param.nkap).^2*(param.B-param.b); % grid more dense at lower end, increase exponents to increase mass at lower end
  
% discretizing AR(1) for income   
    % the process we approximate: log(s_t) = rho*log(s_t-1)+sigmaint*sqrt(1-rho^2)*error_t
    param.nz              = 7;                % number of discretized income states
    param.rho             = 0.2;              % first-order autoregressive coefficient of income
    param.sigmaLR         = 0.4;              % long-run standard deviation of the income process
    param.sigma           = param.sigmaLR*sqrt(1-param.rho^2); % standard deviation of error_t
    
    %Pz is transition matrix of the Markov chain
    %logz is the discretized states of log labor earnings
    %distz is the invariant distribution of Markov chain
    [grid.Pz,grid.logz,grid.distz] = markovappr(param.rho,param.sigma,3,param.nz);
    grid.z               = exp(grid.logz);    % bring back to levels
    param.labor          = grid.z*grid.distz; % aggregate labor is average efficiency units   

    
   % Initialisation
    minrate     = -param.delta;
    maxrate     = (1-param.beta)/(param.beta);  
    err         = 1;
    iter        = 0;
    
% Main Aiyagari loop to iterate on equilibrium interest rate
tic
    while abs(err)>param.tol_r
        iter = iter + 1;
        
    % update guesses
        r0      = 0.5*(maxrate+minrate);                                   % update interest rate guess (bisection)
        k0      = ((r0+param.delta)/(param.alpha*param.A*param.labor^(1-param.alpha)))^(1/(param.alpha-1));    % implied guess for capital stock
        w0      = (1-param.alpha)*(param.A*(param.alpha/(r0+param.delta))^param.alpha)^(1/(1-param.alpha));    % wage implied by interest rate guess
            
    % test whether natural borrowing limit holds
        if (r0>0) && (param.b < -(w0*grid.z(1))/r0)
            fprintf('Warning! Natural borrowing limit violated!\n')
        end
        
        

    
     % solve HH problem for given interest rate
        [cpol,kpol] = solveHH_EGM(r0,w0,0,0,param,grid);
     
     %  compute distribution over households
     %  Note with method 2 is broken, you can fix and email me!
     %  [dist_0] = getDist_method_2(param,grid,kpol,0);
     
     [dist] = getDist_continuous(param,grid,kpol,0);
     % implied aggregate capital supplied by households
        k1 = sum(sum(dist.*kpol));
        r1 = param.alpha*param.A*param.labor^(1-param.alpha)*max(k1,0.01)^(param.alpha-1)-param.delta;
        
     % update interest rate guess
        err = r1-r0;
        if err < 0
            maxrate = r0;
        else
            minrate = r0;
        end
        disp(['k0 = ',num2str(k0),', k1 = ',num2str(k1),', r0 = ',num2str(r0),', r1 = ',num2str(r1)])

    end
  steady.time=toc
  fprintf('Steady State solved, it took me %3.4f seconds. \n\n',[steady.time])
   
% load output structure
 SS.r = r0;
 SS.w = w0;
 SS.k = k0;
 SS.cpol = cpol;
 SS.kpol = kpol;
 SS.dist = dist;
 SS.gridk = grid.k;
 SS.y = param.A*param.labor^(1-param.alpha)*k1^param.alpha;
    
 
 
 
%% Linewidth
   WD=2.5;
%% Fontsize
   FS=14;
    
figure(1)
plot(grid.k,[grid.k' kpol(1,:)' kpol(4,:)' kpol(5,:)'], '-.','linewidth', WD)
ylabel('next period asset holdings','Interpreter','LaTex','Fontsize',FS)
xlabel('current asset holdings','Interpreter','LaTex','Fontsize',FS)
set(gca,'XGrid','off','YGrid','on','Fontsize',FS)
set(gca,'TickLabelInterpreter','LaTex')
xlim([-3 grid.k(param.nkap)])

figure(2)
plot(grid.k,[  sum(dist,1)'],grid.k,[dist(1,:)' dist(3,:)' dist(5,:)'],'-.','linewidth', WD)
ylabel('distribution','Interpreter','LaTex','Fontsize',FS)
xlabel('asset holdings','Interpreter','LaTex','Fontsize',FS)
set(gca,'XGrid','off','YGrid','on','Fontsize',FS)
set(gca,'TickLabelInterpreter','LaTex')
xlim([-3 grid.k(param.nkap)])

figure(3)
plot(grid.k,[cpol(1,:)' cpol(3,:)' cpol(5,:)'], '-.','linewidth', WD)
ylabel('consumption','Interpreter','LaTex','Fontsize',FS)
xlabel('asset holdings','Interpreter','LaTex','Fontsize',FS)
set(gca,'XGrid','off','YGrid','on','Fontsize',FS)
set(gca,'TickLabelInterpreter','LaTex')
xlim([-3 grid.k(param.nkap)])







