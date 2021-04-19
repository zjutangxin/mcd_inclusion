% -------------------------------------------------------------------------
%                      Program Description
% -------------------------------------------------------------------------
%   
% Purpose:
%     - Evaluates all the allocations at a given wage and returns 
%       the excess demand in the capital and labor markets.
%     - Steady state equilibrium.
%     - Buera and Shin (2013, JPE)
%  
% Author:
%     Xin Tang @ International Monetary Fund
%  
% Record of Revisions:
%         Date:                 Description of Changes
%     ============        =================================
%      04/07/2021                  Original Version
% =========================================================================

function varout = fcn_ss_cont(varin)

% ================== Declare global variables =====================
% Shared with main_ss.m
% -------------- Model parameters -----------------
% Technology and preference
global sigma alpha delta beta nu
global eta psi
% Distortions
global tplus tminus q lambda

% -------------- Allocations ----------------------
global polc pola polk poll poly vn ddistss polaind polpi polo
global dpolc dpola dpolk dpoll dpoly dpolo
global egrid edist kgrid kgridc kgrids omega

% -------------- Auxilliary variables ----------------------
global tolv told maxiterv maxiterd
global nkgrid nkgridc sr kint nintk negrid nee tauz kmin kmax
global srsim nsim ademand asupply ldemand lsupply


% =========================================================================
%                      Data dictionary
% =========================================================================
% input variables
w = varin(1) ;
r = varin(2) ;
disp(['r == ', num2str(r), ' w == ', num2str(w)]);
if r+delta <= 0
    r = r+1e-9;
    disp('warning!!!! negative kappa!!!!');
end
    

% computational variables
vnn = zeros(nkgrid,nee) ;
vn = zeros(nkgrid,nee) ;
income = zeros(nkgrid,nee) ;
polk = zeros(nkgrid,nee) ;
poll = zeros(nkgrid,nee) ;
poly = zeros(nkgrid,nee) ;
polpi = zeros(nkgrid,nee) ;
polo = zeros(nkgrid,nee) ;

eps = 10000 ;
loopn = 1 ; 

% =========================================================================
%                   Value function iteration
% =========================================================================
kappa = alpha*w/((1-alpha)*(r+delta)) ;
lstr = (w./(tauz*(1-nu)*(1-alpha)*kappa^(alpha*(1-nu)))).^(-1/nu) ;
kstr = kappa*lstr ;

for indz = 1:1:2
    for inde = 1:1:negrid
        for indk = 1:1:nkgrid
            kcstr = min(kstr(inde,indz),lambda*kgrid(indk)) ;
            lcstr = (w/(tauz(inde,indz)*(1-nu)*(1-alpha)*...
                (kcstr)^(alpha*(1-nu))))^(1/(alpha*(nu-1)-nu)) ;
            pistr = tauz(inde,indz)*(kcstr^alpha*lcstr^(1-alpha))^(1-nu) ...
                - w*lcstr - (r+delta)*kcstr ;
            mop = max(pistr,w) ;     

            polk(indk,(indz-1)*negrid+inde) = kcstr ;
            poll(indk,(indz-1)*negrid+inde) = lcstr ;
            poly(indk,(indz-1)*negrid+inde) = ...
                egrid(inde)*(kcstr^alpha*lcstr^(1-alpha))^(1-nu) ;
            polpi(indk,(indz-1)*negrid+inde) = pistr ;            
            
            if pistr > w
                polo(indk,(indz-1)*negrid+inde) = 1 ;
            end
            income(indk,(indz-1)*negrid+inde) = mop + (1+r)*kgrid(indk) ;
            vn(indk,(indz-1)*negrid+inde) =( income(indk,(indz-1)*negrid+inde)^(1-sigma)-1)/(1-sigma) ;
        end
    end
end

while (eps > tolv) && (loopn < maxiterv)
    if mod(loopn,5) == 0
        disp(['loopn = ', num2str(loopn), '   eps = ', num2str(eps)])
    end
    % initialization
    polaind = nkgridc*ones(nkgrid,nee) ;
    polc = zeros(nkgrid,nee) ;
    pola = zeros(nkgrid,nee) ;
    vtemp = zeros(nkgridc,1) ;
    
    % ================== Occupational Choice =====================
    for indz = 1:1:2
    for inde = 1:1:negrid
        for indk = 1:1:nkgrid
%             kcstr = min(kstr(inde,indz),lambda*kgrid(indk)) ;
%             lcstr = (w/(tauz(inde,indz)*(1-nu)*(1-alpha)*...
%                 (kcstr)^(alpha*(1-nu))))^(1/(alpha*(nu-1)-nu)) ;
%             pistr = tauz(inde,indz)*(kcstr^alpha*lcstr^(1-alpha))^(1-nu) ...
%                 - w*lcstr - (r+delta)*kcstr ;
%             mop = max(pistr,w) ;
%             if pistr > w
%                 polk(indk,(indz-1)*negrid+inde) = kcstr ;
%                 poll(indk,(indz-1)*negrid+inde) = lcstr ;
%                 poly(indk,(indz-1)*negrid+inde) = ...
%                     egrid(inde)*(kcstr^alpha*lcstr^(1-alpha))^(1-nu) ;
%             end
%             income = mop + (1+r)*kgrid(indk) ;
            
            % ================== consumption saving =====================
            incomes = income(indk,(indz-1)*negrid+inde); 
            if incomes <= kint
                nj = max(1,fix((incomes-kmin)*sr*(nintk-1.0)/(kint-kmin)+1.0)-1) ;
            else
                nj = min(nkgridc,fix((nintk-1)*sr+1+...
                    (incomes-kint)/(kmax-kint)*sr*(nkgrid-nintk)))-1 ;
            end
            
            % use monotonicity to accelerate
            nj = min(nj,polaind(min(nkgrid,indk+1),(indz-1)*negrid+inde)) ; 
            for indkc = nj:-1:1
                jj = min(nkgrid-1,fix((indkc-1)/sr)+1) ;
                vtemp(indkc) = psi*...
                    ((vn(jj+1,inde+negrid*(indz-1))-vn(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgridc(indkc)-kgrid(jj))+vn(jj,inde+negrid*(indz-1)));
                vtemp(indkc) = vtemp(indkc) + (1-psi)*...
                    sum(omega'.*((vn(jj+1,:)-vn(jj,:))/(kgrid(jj+1)-kgrid(jj))*(kgridc(indkc)-kgrid(jj))...
                    +vn(jj,:)));
                cons = incomes - kgridc(indkc) ;
                vtemp(indkc) = beta*vtemp(indkc) + (cons^(1-sigma)-1)/(1-sigma) ;
                if (indkc < nj) && (vtemp(indkc) < vtemp(min(indkc+1,nj))) 
                    break ;
                end
                vnn(indk,inde+(indz-1)*negrid) = vtemp(indkc) ;
                pola(indk,inde+(indz-1)*negrid) = kgridc(indkc) ;
                polc(indk,inde+(indz-1)*negrid) = incomes - kgridc(indkc) ;
                polaind(indk,inde+(indz-1)*negrid) = indkc ;
            end
        end % indk
    end % inde
    end % indz
    
    maxiterh = maxiterv ;
    looph = 1;    
    epsh = 1;
    vh = vnn ;
    tolh = tolv ;
    
    % Howard acceleration
    while (looph < maxiterh) && (epsh > tolh)
    for indz = 1:1:2
    for inde = 1:1:negrid
        for indk = 1:1:nkgrid
            aind = polaind(indk,inde+(indz-1)*negrid) ;
            jj = min(nkgrid-1,fix((aind-1)/sr)+1) ;
            vtmp = psi*...
                ((vh(jj+1,inde+negrid*(indz-1))-vh(jj,inde+negrid*(indz-1)))...
                /(kgrid(jj+1)-kgrid(jj))*...
                (kgridc(aind)-kgrid(jj))+vh(jj,inde+negrid*(indz-1)));
            vtmp = vtmp + (1-psi)*...
                sum(omega'.*((vh(jj+1,:)-vh(jj,:))/(kgrid(jj+1)-kgrid(jj))*(kgridc(aind)-kgrid(jj))...
                +vh(jj,:)));
            vnn(indk,inde+negrid*(indz-1)) = (polc(indk,inde+negrid*(indz-1))^(1-sigma)-1)/(1-sigma) ...
                + beta*vtmp ;
        end % indk
    end % inde
    end % indz
    
    epsh = max(max(abs(vnn - vh)));
    vh = vnn ;
    looph = looph + 1;
    end % howard
    
    eps = max(max(abs(vn-vnn)));
    vn = vnn ;
    loopn = loopn + 1 ;
end

% % =========================================================================
% %                   Stationary distribution with interpolation
% % =========================================================================
% % load('./vfi_test_huge1.mat');
eps = 10000;
loopd = 1 ;

% Interpolate policy functions
dpolc = zeros(nsim,nee) ;
dpola = zeros(nsim,nee) ;
dpolk = zeros(nsim,nee) ;
dpoll = zeros(nsim,nee) ;
dpoly = zeros(nsim,nee) ;
dpolo = zeros(nsim,nee) ;
for indz = 1:1:2
for inde = 1:1:negrid
    for indk = 1:1:nsim
        jj = min(nkgrid-1,fix((indk-1)/srsim)+1) ;
        dpolc(indk,inde+(indz-1)*negrid) = ((polc(jj+1,inde+negrid*(indz-1))-polc(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgrids(indk)-kgrid(jj))+polc(jj,inde+negrid*(indz-1)));        
        dpola(indk,inde+(indz-1)*negrid) = ((pola(jj+1,inde+negrid*(indz-1))-pola(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgrids(indk)-kgrid(jj))+pola(jj,inde+negrid*(indz-1)));        
        dpolk(indk,inde+(indz-1)*negrid) = ((polk(jj+1,inde+negrid*(indz-1))-polk(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgrids(indk)-kgrid(jj))+polk(jj,inde+negrid*(indz-1)));        
        dpoll(indk,inde+(indz-1)*negrid) = ((poll(jj+1,inde+negrid*(indz-1))-poll(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgrids(indk)-kgrid(jj))+poll(jj,inde+negrid*(indz-1)));                
        dpoly(indk,inde+(indz-1)*negrid) = ((poly(jj+1,inde+negrid*(indz-1))-poly(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgrids(indk)-kgrid(jj))+poly(jj,inde+negrid*(indz-1)));        
        dpolo(indk,inde+(indz-1)*negrid) = ((polo(jj+1,inde+negrid*(indz-1))-polo(jj,inde+negrid*(indz-1)))...
                    /(kgrid(jj+1)-kgrid(jj))*...
                    (kgrids(indk)-kgrid(jj))+polo(jj,inde+negrid*(indz-1)));                        
    end
end
end

ddistss = 1/nsim*ones(nsim,nee) ;
while (eps > told) && (loopd < maxiterd)
   if mod(loopd,100) == 0
       disp(['loopd == ', num2str(loopd), ' eps == ', num2str(eps)]);
   end
   
   distssn = zeros(nsim,nee) ;
   for indz = 1:1:2
   for inde = 1:1:negrid
       for indk = 1:1:nsim
           if dpola(indk,inde+(indz-1)*negrid) <= kint
               jj = max(1,fix((dpola(indk,inde+(indz-1)*negrid)-kmin)*srsim*(nintk-1.0)/(kint-kmin)+1.0)) ;
           else
               jj = min(nsim-1,fix((nintk-1)*srsim+1+...
                    (dpola(indk,inde+(indz-1)*negrid)-kint)/(kmax-kint)*srsim*(nkgrid-nintk))) ;
           end
%            jj = find(kgrids - dpola(indk,inde+(indz-1)*negrid)>0,1)-1;
           distssn(jj,inde+(indz-1)*negrid) = distssn(jj,inde+(indz-1)*negrid) + ...
               ddistss(indk,inde+(indz-1)*negrid)*(kgrids(jj+1)-dpola(indk,inde+(indz-1)*negrid))/(kgrids(jj+1)-kgrids(jj));
           distssn(jj+1,inde+(indz-1)*negrid) = distssn(jj+1,inde+(indz-1)*negrid) + ...
               ddistss(indk,inde+(indz-1)*negrid)*(dpola(indk,inde+(indz-1)*negrid)-kgrids(jj))/(kgrids(jj+1)-kgrids(jj)) ;
       end
   end
   end
   % marginal distribution of a
   disttemp = ddistss*omega ;
   
   for indj = 1:1:nee
       distssn(:,indj) = psi*distssn(:,indj) + (1-psi)*disttemp ;
   end
   eps = max(max(abs(distssn - ddistss)));
   ddistss = distssn ;
   loopd = loopd + 1;
end

asupply = 0.0 ;
ademand = 0.0 ;
lsupply = 0.0 ;
ldemand = 0.0 ;
dpolo = round(dpolo) ;
maxprof = max(polpi) ;
wgt = zeros(nee,1) ; % worker weight
if maxprof(1) > w
    wgt(1) = w/maxprof(1);
end
if maxprof(negrid+1) > w
    wgt(negrid+1) = w/maxprof(negrid+1);
end
for indz = 1:1:2
    for inde = 1:1:negrid-1
        induse = inde+(indz-1)*negrid;
        if (w > maxprof(induse)) && (w > maxprof(induse+1))
            wgt(induse+1) = 1.0 ;
        elseif (w > maxprof(induse)) && (w < maxprof(induse+1))
            wgt(induse+1) = (w-maxprof(induse))/(maxprof(induse+1)-maxprof(induse)) ;
        end
    end
end

for indz = 1:1:2
for inde = 1:negrid
    if (wgt(inde+(indz-1)*negrid)>0) && (wgt(inde+(indz-1)*negrid) < 1)
        for indk = 1:1:nsim
            if (dpolo(indk,inde+(indz-1)*negrid) > 0)
                ademand = ademand + dpolk(indk,inde+(indz-1)*negrid)*(1-wgt(inde+(indz-1)*negrid))*...
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
                asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*...
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
                ldemand = ldemand + dpoll(indk,inde+(indz-1)*negrid)*(1-wgt(inde+(indz-1)*negrid))*...
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
                lsupply = lsupply + ...
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid)*wgt(inde+(indz-1)*negrid);
            else
                asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*...
                        ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
                lsupply = lsupply + ...
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);                
            end
        end        
    else
        for indk = 1:1:nsim
            if (dpolo(indk,inde+(indz-1)*negrid) > 0)
                ademand = ademand + dpolk(indk,inde+(indz-1)*negrid)*...
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
                asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*...
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);            
                ldemand = ldemand + dpoll(indk,inde+(indz-1)*negrid)*...
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
            else
                asupply = asupply + dpola(indk,inde+(indz-1)*negrid)*...
                    ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
                lsupply = lsupply + ddistss(indk,inde+(indz-1)*negrid)*omega(inde+(indz-1)*negrid);
            end
        end
    end % end if
end
end

% asupply = sum((dpola.*ddistss)*omega) ;
% ademand = sum((dpolk.*ddistss)*omega) ;
% ldemand = sum((dpoll.*ddistss)*omega) ;
% indl = find(dpoll == 0) ;
% dpolsupply = zeros(nsim,nee) ;
% dpolsupply(indl) = 1 ;
% lsupply = sum((dpolsupply.*ddistss)*omega) ;

% intent(out)
varout(1) = ldemand - lsupply ;
varout(2) = ademand - asupply ;

disp('wage and interest ')
disp(varin)
disp('Excess demand: Labor and capital')
disp(varout)

end

% % =========================================================================
% %                   Stationary distribution
% % =========================================================================
% load('./vfi_test_large.mat');
% eps = 10000;
% loopd = 1 ;
% distss = 1/nkgrid*ones(nkgrid,nee) ;
% while (eps > told) && (loopd < maxiterd)
%    if mod(loopd,100) == 0
%        disp(['loopd == ', num2str(loopd)]);
%    end
%    
%    distssn = zeros(nkgrid,nee) ;
%    for indz = 1:1:2
%    for inde = 1:1:negrid
%        for indk = 1:1:nkgrid
%            if pola(indk,inde+(indz-1)*negrid) <= kint
%                jj = fix((pola(indk,inde+(indz-1)*negrid)-kmin)*(nintk-1.)/(kint-kmin)+1.0) ;
%            else
%                jj = min(nkgrid-1,fix(nintk+(pola(indk,inde+(indz-1)*negrid)-kint)/(kmax-kint)*(nkgrid-nintk))) ;
%            end
%            distssn(jj,inde+(indz-1)*negrid) = distssn(jj,inde+(indz-1)*negrid) + ...
%                distss(indk,inde+(indz-1)*negrid)*(kgrid(jj+1)-pola(indk,inde+(indz-1)*negrid))/(kgrid(jj+1)-kgrid(jj));
%            distssn(jj+1,inde+(indz-1)*negrid) = distssn(jj+1,inde+(indz-1)*negrid) + ...
%                distss(indk,inde+(indz-1)*negrid)*(pola(indk,inde+(indz-1)*negrid)-kgrid(jj))/(kgrid(jj+1)-kgrid(jj)) ;
%        end
%    end
%    end
%    % marginal distribution of a
%    disttemp = distss*omega ;
%    
%    for indj = 1:1:nee
%        distssn(:,indj) = psi*distssn(:,indj) + (1-psi)*disttemp ;
%    end
%    eps = max(max(abs(distssn - distss)));
%    distss = distssn ;
% end

% =========================================================================
%                   Market excess demand
% =========================================================================
% asupply = sum((pola.*distss)*omega) ;
% ademand = sum((polk.*distss)*omega) ;
% ldemand = sum((poll.*distss)*omega) ;
% indl = find(poll == 0) ;
% polsupply = zeros(nkgrid,nee) ;
% polsupply(indl) = 1 ;
% lsupply = sum((polsupply.*distss)*omega) ;
