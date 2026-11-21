
% RENEWAL_BIFURCATION  Three-pathway CTMC: scalar steady-state diagram.
%
% States: Q -> F -> C_1 -> ... -> C_n -> Q.
% Extra edges: C_k -> F and C_l -> Y_1 -> ... -> Y_m -> Q.
% The k and l indices can be in either order, or coincide.
%
% Parameterize equilibria by z:
%   p_F = H(z) = 1/(gamma0*meanCycleTime(z)), u = z - a*H(z).
% Scalar continuous-time slope criterion:
%   stable   if a*H'(z) < 1, equivalently du/dz > 0;
%   unstable if a*H'(z) > 1, equivalently du/dz < 0.
% This is the criterion for dp_F/dt = H(a*p_F + u) - p_F,
% not the criterion for discrete fixed-point iteration.
% No toolboxes required. Parameters below are illustrative, not fitted.

%% Parameters
% Canonical pathway
params.n       = 10;
params.gamma0  = 1;
params.gamma   = 1*ones(1,params.n);
params.Amax    = 1;
params.muA     = 3.5;
params.thetaA  = 2.5;
params.a       = 8;

% Recurrence pathway C_k -> F, with rate beta_k P(z)
params.k       = 2;
params.beta    = 8;
params.muP     = 2.5;
params.thetaP  = 3.5;

% Sequestration pathway 
params.kappa = 0.1;
params.l       = 1;
params.m       = 30;
params.alpha   = 1*ones(1,params.m);  % Y_j forward rates, including Y_m -> Q_S

zRange = [-2, 12];               % Input z*, not external input u
nGrid = 8001;                    % Increase for narrow branch features

%% Additional panels
uSelected = 1.7;                % External input at which to show distributions
tMax = 100;                      % [] chooses a shared time horizon for all panels
nTime = 4001;                   % Number of time samples per fixed point
survivalFloor = 1e-10;          % Omit hazards once total survival is negligible

% Axis scales for additional panels only: 'linear' or 'log'.
timeXScale = 'linear';          % Shared by every density and hazard panel
densityYScale = 'linear';
hazardYScale = 'linear';
logYDecades = 8;                % Log y-limits extend this far below group peak

canonicalColor      = [0.3,0.8,0.3];
recurrenceColor     = [1.00,0.30,0.70];
sequestrationColor  = [0.00,0.4,0.70];
totalColor          = [0,0,0];
componentLineWidth  = 1.5;
totalLineWidth      = 2;

%% Check parameters
assert(params.n >= 1 && params.n == fix(params.n));
assert(params.m >= 1 && params.m == fix(params.m));
assert(all([params.k,params.l] >= 1 & [params.k,params.l] <= params.n));
assert(all([params.k,params.l] == fix([params.k,params.l])));
assert(numel(params.gamma) == params.n && all(params.gamma > 0));
assert(numel(params.alpha) == params.m && all(params.alpha > 0));
assert(params.gamma0 > 0 && params.Amax > 0);
assert(params.beta >= 0 && params.kappa >= 0);
assert(all([params.thetaA,params.thetaP] > 0));
assert(zRange(1) < zRange(2) && nGrid >= 3);
assert(isscalar(uSelected) && isfinite(uSelected));
assert(isempty(tMax) || (isscalar(tMax) && isfinite(tMax) && tMax > 0));
assert(nTime >= 3 && nTime == fix(nTime));
assert(survivalFloor > 0 && survivalFloor < 1);
assert(any(strcmp(timeXScale,{'linear','log'})));
assert(any(strcmp(densityYScale,{'linear','log'})));
assert(any(strcmp(hazardYScale,{'linear','log'})));
assert(isscalar(logYDecades) && isfinite(logYDecades) && logYDecades > 0);

%% Parameterize the entire equilibrium curve, including folds
z = linspace(zRange(1),zRange(2),nGrid);
[pF,~,~] = mean_renewal_response(z,params);
Hp = gradient(pF,z);
duDz = 1 - params.a*Hp;

% Locate turning points and insert them as shared branch endpoints.
% A dense z grid resolves branches even when u(z) is not monotone.
crossings = find(duDz(1:end-1).*duDz(2:end) < 0);
zFold = z(duDz == 0);
for j = 1:numel(crossings)
    i = crossings(j);
    zFold(end+1) = fzero(@(x) slope_u(x,params),[z(i),z(i+1)]); 
end
zFold = unique(zFold);
z = unique([z,zFold]);            % Keep z order; DO NOT sort by u

[pF,mu_array,rho_array] = mean_renewal_response(z,params);
Hp = gradient(pF,z);
meanCycleTime = mu_array(:,4).';
rho = struct('C',rho_array(:,1).','R',rho_array(:,2).','S',rho_array(:,3).');
u = z - params.a*pF;
duDz = 1 - params.a*Hp;
isFold = ismember(z,zFold);
isStable = duDz > 0 & ~isFold;
isUnstable = duDz < 0 & ~isFold;

% NaNs prevent lines joining disconnected branch segments.
uStable = u;    yStable = pF;
uUnstable = u;  yUnstable = pF;
uStable(~(isStable | isFold)) = NaN;
yStable(~(isStable | isFold)) = NaN;
uUnstable(~(isUnstable | isFold)) = NaN;
yUnstable(~(isUnstable | isFold)) = NaN;

%% Stable fixed points at the selected external input
residual = u - uSelected;
rootIntervals = find(residual(1:end-1).*residual(2:end) < 0);
zSelected = z(abs(residual) < 1e-10);
for j = 1:numel(rootIntervals)
    i = rootIntervals(j);
    zSelected(end+1) = fzero(@(x) x ...
        - params.a*mean_renewal_response(x,params) - uSelected, ...
        [z(i),z(i+1)]);
end
zSelected = sort(zSelected);
if ~isempty(zSelected)
    zSelected = zSelected([true,diff(zSelected) > 1e-8]);
end
selectedSlope = arrayfun(@(x) slope_u(x,params),zSelected);
zSelected = zSelected(selectedSlope > 1e-7); % Exclude unstable/neutral roots
[pSelected,~,~] = mean_renewal_response(zSelected,params);
nStable = numel(zSelected);
returnTime = cell(nStable,1);
pointLabels = repmat({'Stable state'},nStable,1);
if nStable == 2
    [~,activityOrder] = sort(pSelected);
    pointLabels{activityOrder(1)} = 'Low-activity';
    pointLabels{activityOrder(2)} = 'High-activity';
end

% Precompute path rates so every lower panel covers the same time range.
selectedRates = cell(nStable,1);
pathHorizons = zeros(nStable,3);
for j = 1:nStable
    A = params.Amax*sigmoid(zSelected(j),params.muA,params.thetaA);
    P = sigmoid(zSelected(j),params.muP,params.thetaP);
    lambda = params.gamma;
    lambda(params.k) = lambda(params.k) + params.beta*P;
    lambda(params.l) = lambda(params.l) + params.kappa;
    selectedRates{j} = {[params.gamma0,lambda,A], ...
        [params.gamma0,lambda(1:params.k)], ...
        [params.gamma0,lambda(1:params.l),params.alpha,A]};
    pathHorizons(j,:) = cellfun(@(v) sum(1./v) ...
        + 12*sqrt(sum(1./v.^2)),selectedRates{j});
end
if isempty(tMax)
    sharedTimeMax = max(pathHorizons(:));
else
    sharedTimeMax = tMax;
end
densityAxes = gobjects(nStable,1);
hazardAxes = gobjects(nStable,1);

%% Panel layout: bifurcation above, density and hazard below
figure('Color','w','Position',[100,100,900,240*(1+nStable)]);
panelLayout = tiledlayout(1+nStable,2,'TileSpacing','compact','Padding','compact');
nexttile(panelLayout,[1,2]);
plot(uStable,yStable,'r-','LineWidth',2); hold on;
plot(uUnstable,yUnstable,'k:','LineWidth',2);
xline(uSelected,'--','Color',[0,0,1],'LineWidth',1.5,'HandleVisibility','off');
scatter(uSelected*ones(size(pSelected)),pSelected,45,'o','MarkerFaceColor','w', ...
    'MarkerEdgeColor',[0 0 1],'LineWidth',1,'HandleVisibility','off');
xlabel('$u = z^* - a p_F^*$');
ylabel('$p_F^*$');
title('Three-pathway renewal model: steady-state bifurcation diagram');
legend('Stable','Unstable','Location','best');
grid on; box on;
set(gca,'FontSize',12,'XLim',[0 6],'YLim',[0 1.2*max(yStable,[],"all")]);

if nStable == 0
    warning('No strictly stable fixed point found at u = %g in zRange.',uSelected);
end

for j = 1:nStable
    % Conditional return times are sums of exponential holding times.
    rates = selectedRates{j};
    [~,~,weights] = mean_renewal_response(zSelected(j),params);
    rhoText = sprintf('%s %.3g, %s %.3g, %s %.3g', ...
        '$\rho^{C}=$',weights(1),'$\rho^{R}=$',weights(2), ...
        '$\rho^{S}=$',weights(3));

    % Resolve both short recurrence times and long sequestration times.
    pathHorizon = pathHorizons(j,:);
    tEnd = sharedTimeMax;
    tBreak = min(pathHorizon);
    if tBreak < tEnd/2
        t = unique([linspace(0,tBreak,ceil(nTime/2)), ...
            linspace(tBreak,tEnd,floor(nTime/2)+1)]).';
    else
        t = linspace(0,tEnd,nTime).';
    end

    conditionalPDF = zeros(numel(t),3);
    conditionalSurvival = zeros(numel(t),3);
    for pathway = 1:3
        [conditionalPDF(:,pathway),conditionalSurvival(:,pathway)] = ...
            phase_type_curve(rates{pathway},t);
    end
    densityComponents = conditionalPDF.*weights;
    survivalComponents = conditionalSurvival.*weights;
    totalPDF = sum(densityComponents,2);
    totalSurvival = sum(survivalComponents,2);
    hazardComponents = NaN(size(densityComponents));
    valid = isfinite(totalSurvival) & totalSurvival > 0;
    % rho_j*f_j/S = (rho_j*S_j/S)*(f_j/S_j): survival-weighted hazards.
    hazardComponents(valid,:) = densityComponents(valid,:)./totalSurvival(valid);
    totalHazard = sum(hazardComponents,2);

    densityAxes(j) = nexttile(panelLayout,2*j+1);
    plot(t,densityComponents(:,1),'-','Color',canonicalColor, ...
        'LineWidth',componentLineWidth); hold on;
    plot(t,densityComponents(:,2),'-','Color',recurrenceColor, ...
        'LineWidth',componentLineWidth);
    plot(t,densityComponents(:,3),'-','Color',sequestrationColor, ...
        'LineWidth',componentLineWidth);
    plot(t,totalPDF,':','Color',totalColor,'LineWidth',totalLineWidth);
    xlabel('Return time $t$'); ylabel('Return-time density $f(t)$');
    title(sprintf('%s: $u = %.4g$, $p_F^* = %.4g$, $z^* = %.4g$', ...
        pointLabels{j},uSelected,pSelected(j),zSelected(j)));
    text(0.98,0.95,rhoText,'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'BackgroundColor','none','Margin',2,'FontSize',11,'Interpreter','latex');
    grid on; box on; set(gca,'FontSize',12); xlim([0,tEnd]);

    hazardAxes(j) = nexttile(panelLayout,2*j+2);
    plot(t,hazardComponents(:,1),'-','Color',canonicalColor, ...
        'LineWidth',componentLineWidth); hold on;
    plot(t,hazardComponents(:,2),'-','Color',recurrenceColor, ...
        'LineWidth',componentLineWidth);
    plot(t,hazardComponents(:,3),'-','Color',sequestrationColor, ...
        'LineWidth',componentLineWidth);
    plot(t,totalHazard,':','Color',totalColor,'LineWidth',totalLineWidth);
    xlabel('Return time $t$'); ylabel('Hazard $h(t)$');
    title(sprintf('%s: survival-weighted hazard, $p_F^* = %.4g$', ...
        pointLabels{j},pSelected(j)));
    text(0.98,0.95,rhoText,'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'BackgroundColor','none','Margin',2,'FontSize',11,'Interpreter','latex');
    grid on; box on; set(gca,'FontSize',12); xlim([0,tEnd]);

    returnTime{j} = struct('z',zSelected(j),'pF',pSelected(j),'t',t, ...
        'weights',weights,'conditionalPDF',conditionalPDF, ...
        'conditionalSurvival',conditionalSurvival, ...
        'densityComponents',densityComponents,'totalPDF',totalPDF, ...
        'survivalComponents',survivalComponents,'totalSurvival',totalSurvival, ...
        'hazardComponents',hazardComponents,'totalHazard',totalHazard);
end

% Shared time limits, and separate shared density and hazard ordinate limits.
if nStable > 0
    densityMax = max(cellfun(@(r) max(r.totalPDF),returnTime));
    hazardMax = max(cellfun(@(r) ...
        max(r.totalHazard(isfinite(r.totalHazard))),returnTime));
    densityPeak = max(densityMax,eps);
    hazardPeak = max(hazardMax,eps);
    densityLower = 0;
    hazardLower = 0;
    if strcmp(densityYScale,'log')
        densityLower = max(realmin,densityPeak*10^(-logYDecades));
    end
    if strcmp(hazardYScale,'log')
        hazardLower = max(realmin,hazardPeak*10^(-logYDecades));
    end
    set(densityAxes,'YLim',[densityLower,1.05*densityPeak], ...
        'YScale',densityYScale);
    set(hazardAxes,'YLim',[hazardLower,1.05*hazardPeak], ...
        'YScale',hazardYScale);
    timeAxes = [densityAxes;hazardAxes];
    timeLower = 0;
    if strcmp(timeXScale,'log')
        timeLower = min(cellfun(@(r) min(r.t(r.t > 0)),returnTime));
    end
    set(timeAxes,'XLim',[timeLower,sharedTimeMax],'XScale',timeXScale);
    linkaxes(timeAxes,'x');

    % Mark each mean return time at the density panel's lower boundary.
    for j = 1:nStable
        meanReturnTime = 1/(params.gamma0*pSelected(j));
        if meanReturnTime >= timeLower && meanReturnTime <= sharedTimeMax
            plot(densityAxes(j),meanReturnTime,densityLower,'*', ...
                'Color',[0.8 0 0],'MarkerSize',10,'LineWidth',1.5, ...
                'Clipping','off','HandleVisibility','off');
        end
    end

    % One shared pathway legend for all additional panels.
    componentLegend = legend(densityAxes(1), ...
        'Canonical','Recurrence','Sequestration','Total', ...
        'Orientation','horizontal','NumColumns',4);
    componentLegend.Layout.Tile = 'south';
    componentLegend.FontSize = 11;
end

% Apply LaTeX to titles, labels, legends, and numeric tick labels.
set(findall(gcf,'Type','text'),'Interpreter','latex');
set(findall(gcf,'Type','legend'),'Interpreter','latex');
set(findall(gcf,'Type','axes'),'TickLabelInterpreter','latex');

%% Return computed values for further analysis
out.parameters = params;
out.z = z;
out.u = u;
out.pF = pF;
out.meanCycleTime = meanCycleTime;
out.dpFdz = Hp;
out.duDz = duDz;
out.stable = isStable;
out.unstable = isUnstable;
out.fold = isFold;               % Neutral points, included in both lines
out.pathwayProbabilities = rho;  % Fields C, R, S
out.turningPoints = table(z(isFold).',u(isFold).',pF(isFold).', ...
    'VariableNames',{'z','u','pF'});
disp(out.turningPoints);
out.selectedInput = uSelected;
out.selectedStableStates = table(zSelected(:),pSelected(:), ...
    'VariableNames',{'z','pF'});
out.returnTime = returnTime;


function [pF,mu_array,rho_array] = mean_renewal_response(z,params)
% Exact renewal mean at prescribed constant z.
% At C_i, all outgoing hazards add, including when k == l.

z = z(:);
A = sigmoid(z,params.muA,params.thetaA);
P = sigmoid(z,params.muP,params.thetaP);

lambda = repmat(params.gamma,numel(z),1);
lambda(:,params.k) = lambda(:,params.k) + params.beta*P;
lambda(:,params.l) = lambda(:,params.l) + params.kappa;
disp(size(lambda));

% Mean return time for each pathway
mu_C = 1/params.gamma0 + sum(1./lambda,2) + 1./(params.Amax*A);
mu_R = 1/params.gamma0 + sum(1./lambda(:,1:params.k),2);
mu_S = 1/params.gamma0 + sum(1./lambda(:,1:params.l),2) + sum(1./params.alpha) + 1./(params.Amax*A);

% Pathway probabilites
canonical = params.gamma ./ lambda;
reach = cumprod([ones(numel(z),1),canonical(:,1:end-1)],2);
rho_R = reach(:,params.k).*params.beta.*P ./ lambda(:,params.k);
rho_S = reach(:,params.l).*params.kappa ./ lambda(:,params.l);
rho_C = 1 - rho_R - rho_S;

% Total mean return time
mu = rho_C.*mu_C + rho_R.*mu_R + rho_S.*mu_S;

% Firing occupancy
pF = 1./(params.gamma0*mu);
pF = pF.';                     % Match the row vector z in the script

% Output
mu_array = [mu_C,mu_R,mu_S,mu];
rho_array = [rho_C,rho_R,rho_S];

end

function value = slope_u(z,params)
h = 1e-5*max(1,abs(z));
[pF,~,~] = mean_renewal_response([z-h,z+h],params);
value = 1 - params.a*(pF(2)-pF(1))/(2*h);
end

function s = sigmoid(z,mu,theta)

s = 1 + exp(-theta*(z-mu));
s = 1./s;

end

function [density,survival] = phase_type_curve(rates,t)
% Exact phase-type formula: e_1*expm(T*t), for a serial exponential chain.
% Matrix propagation avoids partial-fraction singularities at repeated rates.
rates = rates(:);
nStages = numel(rates);
T = diag(-rates) + diag(rates(1:end-1),1);
state = zeros(1,nStages);
state(1) = 1;
density = zeros(numel(t),1);
survival = zeros(numel(t),1);
density(1) = state(end)*rates(end);
survival(1) = 1;
lastStep = NaN;
for i = 2:numel(t)
    dt = t(i)-t(i-1);
    % The time grid has at most two uniform segments: reuse propagators.
    if isnan(lastStep) || abs(dt-lastStep) > 1e-9*max(dt,lastStep)
        propagator = expm(T*dt);
        lastStep = dt;
    end
    state = state*propagator;
    density(i) = state(end)*rates(end);
    survival(i) = sum(state);
end
end
