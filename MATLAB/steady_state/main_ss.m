% -------------------------------------------------------------------------
%                      Program Description
% -------------------------------------------------------------------------
%   
% Purpose:
%     - Main program to solve the stationary equilibrium. 
%     - Steady state.
%     - Buera and Shin (2013, JPE)
%  
% Author:
%     - Xin Tang @ International Monetary Fund
% 
% File Dependence:
%     - fcn.m: evaluates all the allocations at a given wage and returns 
%       the excess demand in the capital and labor markets.
%  
% Record of Revisions:
%         Date:                 Description of Changes
%     ============        =================================
%      04/07/2021                 Original Version
% =========================================================================

clc ;
clear all; 
tic;
% ================== Declare global variables =====================
% Shared with fcn.m
% -------------- Model parameters -----------------
% Technology and preference
global sigma alpha delta beta nu
global eta psi
% Distortions
global tplus tminus q lambda

% -------------- Allocations ----------------------
global polc pola polk poll poly vn ddistss polaind polpi polo
global dpolc dpola dpolk dpoll dpoly dpolo dpolpi
global egrid edist kgrid kgridc kgrids omega wgt maxprof

% -------------- Auxilliary variables ----------------------
global tolv told maxiterv maxiterd
global nkgrid nkgridc sr kint nintk negrid nee tauz kmin kmax
global srsim nsim ldemand lsupply ademand asupply

% =========================================================================
%                      Data dictionary
% =========================================================================
% Utility variables
tolv = 1e-6 ;
told = 1e-7 ;
tolfsolve = 1e-2 ;
maxiterv = 1000 ;
maxiterd = 1000 ;

% Technology and preference
beta  = 0.9040 ;    % discount factor
sigma = 1.5 ;       % CRRA
alpha = 1.0/3 ;     % capital share
delta = 0.06 ;      % depreciation rate
nu    = 0.21 ;      % span-of-control
eta   = 4.15 ;      % Pareto tail parameter
psi   = 0.8940 ;    % persistence of e-tau

% Distortions
lambda = 1.35 ;     % financial frictions
tplus  = 0.50 ;     % tax
tminus = -0.2975 ;  % subsidy
q      = 7.0 ;      % tax probability Pr = 1-exp(-qe)
% TODO: need to verify tplus and tminus

% Computational parameters
kmax    = 12500 ;
kmin    = 1e-9 ;
nkgrid  = 501 ;
sr      = 10 ;
nkgridc = sr*(nkgrid-1)+1 ;
srsim = 15 ;
nsim = srsim*(nkgrid-1)+1 ;
kint    = 25.0 ;
nintk   = 201 ;
negrid  = 12 ;
nee = negrid*2 ;

% Computational variables
kgrid  = zeros(nkgrid,1) ;
kgridc = zeros(nkgridc,1) ;
kgrids = zeros(nsim,1) ;
egrid  = zeros(negrid,1) ;
edist  = zeros(negrid,1) ;
omega  = zeros(nee,1) ;   % joint pdf of e-tau
tauz = zeros(negrid,2) ;  % (1-tau)*e

% Allocations
polc = zeros(nkgrid,nee) ;
pola = zeros(nkgrid,nee) ;  % savings
polaind = zeros(nkgrid,nee) ;  % savings index
polk = zeros(nkgrid,nee) ;  % capital demand
poll = zeros(nkgrid,nee) ;
poly = zeros(nkgrid,nee) ;
polo = zeros(nkgrid,nee) ;
vn   = zeros(nkgrid,nee) ;  % value function
ddistss = zeros(nsim,nee) ;

% =========================================================================
%                   Discretization of state space
% =========================================================================
% Capital grid
for indi = 1:1:nintk
    kgrid(indi) = kmin+(indi-1.0)*(kint-kmin)/(nintk-1.0);
end

for indi = nintk+1:1:nkgrid
    kgrid(indi) = kint+(indi-nintk)*(kmax-kint)/(nkgrid-nintk);
end

for indi = 1:1:nkgrid-1
    for indj = 1:1:sr+1
        kgridc((indi-1)*sr+indj) = kgrid(indi)+...
            (indj-1)*(kgrid(indi+1)-kgrid(indi))/sr;
    end
end

for indi = 1:1:nkgrid-1
    for indj = 1:1:srsim+1
        kgrids((indi-1)*srsim+indj) = kgrid(indi)+...
            (indj-1)*(kgrid(indi+1)-kgrid(indi))/srsim;
    end
end

% Entrepreneur productivity
egrid(negrid)   = exp(-log(1-0.9995)/eta) ;
egrid(negrid-1) = exp(-log(1-0.999)/eta) ;
egrid(negrid-2) = exp(-log(1-0.998)/eta) ;
egrid(1)        = exp(-log(1-0.3670)/eta) ;
incr = (egrid(negrid-2)-egrid(1))/(negrid-2-1) ;
for inde = 2:1:negrid-2
    egrid(inde) = egrid(inde-1)+incr ;
end
me = 1-egrid.^-eta ;
for inde = 2:1:negrid
    edist(inde) = (me(inde) - me(inde-1))/me(negrid) ;
end
edist(1) = me(1)/me(negrid) ;
% verify mean and variance
% mean = eta/(eta-1), var = (1/(eta-1))^2*(eta/(eta-2)) 
emean = eta/(eta-1);
evar = (1/(eta-1))^2*(eta/(eta-2)) ;
disp('Pareto Distribution: Compare Mean and Variance')
disp('-------------------------')
disp('Theoretical: mean and var')
disp([emean evar])
disp('Numerical: mean and var')
nummean = sum(egrid.*edist) ;
numvar = sum((egrid-nummean).^2.*edist) ;
disp([nummean numvar])
egrid = 0.2*egrid ;     % Use log as the TFP

omega(1:negrid,1) = edist.*(1-exp(-q*egrid)) ;
omega(negrid+1:nee,1) = edist.*(exp(-q*egrid)) ;
tauz(:,1) = (1-tplus)*egrid ;
tauz(:,2) = (1-tminus)*egrid ;

% =========================================================================
%                   Solve Equilibrium Prices
% =========================================================================
% initial guess
% wguess = 0.0803852580114;
% rguess = -0.050125 ;
wguess = 8.5284241397084265E-002 ;
rguess = -4.8792179372680620E-002 ;
varin = [wguess rguess] ;

% evaluate excess demand
% fval = fcn_ss(varin) ;
% fval = fcn_ss_cont(varin) ;

% % bisection
% wbot = 0.0853852580114 ;
% wup = 0.0873852580114 ;
% wguess = (wbot+wup)/2 ;
% diff = 1 ;
% loopw = 1 ;
% while (diff > 0.01) && (loopw <= 20)
% %    fval = fcn_ss([wguess rguess]) ;
%    fval = fcn_ss_cont([wguess rguess]) ;
%    diff = abs(fval(1)) ;
%    if fval(1) > 0
%        wbot = wguess;
%    else
%        wup = wguess;
%    end
%    wguess = (wbot+wup)/2 ;
%    loopw = loopw + 1;
% end

% solve for equilibrium price
% [eprice fval] = fsolve(@fcn_ss_cont,[wguess rguess],...
%     optimoptions('fsolve','Display','iter',...
%     'OptimalityTolerance',tolfsolve,'MaxIterations',4000,...
%     'FiniteDifferenceStepSize',1e-3));
[eprice fval] = fsolve(@fcn_ss_cont,[wguess rguess],...
    optimoptions('fsolve','Display','iter',...
    'OptimalityTolerance',tolfsolve,'MaxIterations',4000,...
    'FiniteDifferenceStepSize',1e-3));
wss = eprice(1) ;
rss = eprice(2) ;

time = toc ;
save ss_hybrd1.mat;
% save ss_eval.mat ;

disp('total running')
disp(time) ;