%% Canonical CTMC: bifurcation, return-time density, and hazard
% The canonical cycle is
%   F --gamma0--> C1 --gamma--> ... --gamma--> Cn --gamma--> Q
%     --A(z)--> F,
% with z = w*pF + u and total occupancy equal to one.
%
% The return time T_FF starts on entry to F and ends on the next entry to F.
% Its distribution is evaluated exactly as a phase-type distribution.

clear; close all; clc;
save_option = 1;

%% User controls
base.n       = 10;      % Number of recovery states C_1,...,C_n
base.gamma0  = 3;       % F -> C_1 rate (set separately)
base.gamma   = 3;       % Common rate for all remaining canonical gamma_i
base.Amax    = 1;
base.muA     = 3.5;
base.thetaA  = 3.5;
base.w       = 8;
base.nu      = 1;

% Row-specific drives for the return-time and hazard curves. 
uEvalByRow        = [3 2.7 3 2.4 3.25 2.5];
equilibriumBranch = 'high';  % 'high' or 'low'

uLimits = [0, 4];
zGrid   = linspace(-4, 14, 30001);
nTime   = 500;

% Independent x-axis scales for the return-time density and hazard panels.
returnTimeXScale = "linear";  % "linear" or "log"
hazardXScale     = "linear";  % "linear" or "log"

% The calculation grid supports both display modes.
timeEvalLimits  = [0,20];
timePositiveMin = 1e-3;

% Limits and ticks selected automatically according to each scale.
linearTimeXLim        = [0,20];
linearTimeXTicks      = 0:5:20;
linearTimeXTickLabels = ["0","5","10","15","20"];

logTimeXLim        = [0.01,2];
logTimeXTicks      = [0.01,0.1,1];
logTimeXTickLabels = ["0.01","0.1","1"];

returnTimeAxis = makeTimeAxisSettings(returnTimeXScale, ...
    linearTimeXLim,linearTimeXTicks,linearTimeXTickLabels, ...
    logTimeXLim,logTimeXTicks,logTimeXTickLabels);
hazardAxis = makeTimeAxisSettings(hazardXScale, ...
    linearTimeXLim,linearTimeXTicks,linearTimeXTickLabels, ...
    logTimeXLim,logTimeXTicks,logTimeXTickLabels);

% Parameter sweeps. Edit these values as required.
sweeps(1) = makeSweep('gamma',  [2 4 8 16],       '$\gamma$');
sweeps(2) = makeSweep('Amax',   [0.5 1 2 4],       '$\mathcal{A}_{\max}$');
sweeps(3) = makeSweep('w',      [3, 6, 9, 12],        '$a$');
sweeps(4) = makeSweep('muA',    [4, 3.5, 3, 2.5], '$\mu_A$');
sweeps(5) = makeSweep('thetaA', [1.5, 2.5, 4.5, 6.5], '$\theta_A$');
sweeps(6) = makeSweep('nu',     [1, 1.5, 2, 3],       '$\nu$');

colors = [0.5 0.2 1; 0 0.5 1; 0 0.5 0; 0.8 0.4 0];

if numel(uEvalByRow) ~= numel(sweeps)
    error('uEvalByRow must contain one value for each sweep row.');
end
if timeEvalLimits(1) < 0 || timeEvalLimits(2) <= timeEvalLimits(1)
    error('timeEvalLimits must be nonnegative and strictly increasing.');
end
if timePositiveMin <= 0 || timePositiveMin >= timeEvalLimits(2)
    error('timePositiveMin must lie strictly between zero and timeEvalLimits(2).');
end
validateTimeAxisSettings(returnTimeAxis,timeEvalLimits,"returnTimeAxis");
validateTimeAxisSettings(hazardAxis,timeEvalLimits,"hazardAxis");

%% Precompute every row and construct one shared logarithmic time grid
curveRows = cell(1,numel(sweeps));
for row = 1:numel(sweeps)
    sweep = sweeps(row);
    uEval = uEvalByRow(row);
    nCurves = numel(sweep.values);
    curve = repmat(struct(),1,nCurves);
    meanReturn = zeros(1,nCurves);

    for j = 1:nCurves
        p = base;
        p.(sweep.field) = sweep.values(j);

        curve(j).p = p;
        symbolCore = sweep.symbol(2:end-1);
        curve(j).label = num2str(sweep.values(j));
        curve(j).vary = "$"+symbolCore+"$";
        curve(j).bif = bifurcationData(p,zGrid);
        curve(j).eq = selectStableEquilibrium( ...
            p,uEval,zGrid,equilibriumBranch);
        meanReturn(j) = 1/p.gamma0 + p.n/p.gamma + 1/curve(j).eq.A;
        curve(j).meanReturn = meanReturn(j);
    end

    curveRows{row} = curve;
end

t = [0;logspace(log10(timePositiveMin), ...
    log10(timeEvalLimits(2)),nTime-1).'];

% Evaluate all time-dependent curves before plotting so every density panel
% can use exactly the same vertical limits.
globalPdfMax = 0;
for row = 1:numel(sweeps)
    curve = curveRows{row};
    for j = 1:numel(curve)
        [curve(j).pdf,curve(j).hazard] = returnTimeStatistics( ...
            curve(j).p,curve(j).eq.A,t);
        globalPdfMax = max(globalPdfMax,max(curve(j).pdf));
    end
    curveRows{row} = curve;
end
pdfLimits = [0,1.06*globalPdfMax];

%% Figure
fig = figure('Color','w','Units','pixels','Position',[40, 40, 900, 1600]);
tl = tiledlayout(fig,6,4,'TileSpacing','compact','Padding','compact');

returnTimeAxes = gobjects(1,numel(sweeps));
hazardAxes = gobjects(1,numel(sweeps));
for row = 1:numel(sweeps)
    sweep = sweeps(row);
    uEval = uEvalByRow(row);
    curve = curveRows{row};
    nCurves = numel(curve);

    % Columns 1-2: one-parameter bifurcation diagrams.
    ax1 = nexttile(tl,[1,2]);
    plotBifurcationPanel(ax1,curve,colors,uLimits,uEval);
    addRowLabel(ax1,row);

    % Column 3: first-return-time probability density.
    ax2 = nexttile(tl);
    hold(ax2,'on');
    for j = 1:nCurves
        plot(ax2,t,curve(j).pdf,'Color',colors(j,:), ...
            'LineWidth',1.6,'HandleVisibility','off');
    end
    formatAxis(ax2);
    xlabel(ax2,'$t$','Interpreter','latex');
    ylabel(ax2,'$f(t)$','Interpreter','latex');
    set(ax2,"XScale",returnTimeAxis.xScale, ...
        "XLim",returnTimeAxis.xLim, ...
        "XTick",returnTimeAxis.xTicks, ...
        "XTickLabel",returnTimeAxis.xTickLabels);
    set(ax2,'YLim',[0.001,2],'YScale','log','YTick',[0.001 0.01 0.1 1]);
    timeLimits = xlim(ax2);
    densityLimits = ylim(ax2);
    for j = 1:nCurves
        meanTime = curve(j).meanReturn;
        if meanTime >= timeLimits(1) && meanTime <= timeLimits(2)
            plot(ax2,meanTime,densityLimits(1),'*','Color',colors(j,:), ...
                'MarkerSize',10,'LineWidth',1.5,'Clipping','off', ...
                'HandleVisibility','off');
        end
    end
    returnTimeAxes(row) = ax2;

    % Column 4: hazard of the first-return-time distribution.
    ax3 = nexttile(tl);
    hold(ax3,'on');
    for j = 1:nCurves
        plot(ax3,t,curve(j).hazard,'Color',colors(j,:), ...
            'LineWidth',1.6,'HandleVisibility','off');
    end
    formatAxis(ax3);
    xlabel(ax3,'$t$','Interpreter','latex');
    ylabel(ax3,'$h(t)$','Interpreter','latex');
    set(ax3,"XScale",hazardAxis.xScale, ...
        "XLim",hazardAxis.xLim, ...
        "XTick",hazardAxis.xTicks, ...
        "XTickLabel",hazardAxis.xTickLabels);
    set(ax3,'YLim',[0.01,2],'YScale','log','YTick',[0.01 0.1 1]);
    hazardAxes(row) = ax3;

    % Column titles appear only on the top row.
    if row == 1
        title(ax1,'Bifurcation diagram', ...
            'Interpreter','latex','FontSize',12);
        title(ax2,'First return time distribution', ...
            'Interpreter','latex','FontSize',11);
        title(ax3,'Hazard function', ...
            'Interpreter','latex','FontSize',12);
    end
end

% Link each column independently so the two x-scales can differ.
linkaxes(returnTimeAxes,"x");
linkaxes(hazardAxes,"x");
set(returnTimeAxes,"XScale",returnTimeAxis.xScale, ...
    "XLim",returnTimeAxis.xLim,"XTick",returnTimeAxis.xTicks, ...
    "XTickLabel",returnTimeAxis.xTickLabels);
set(hazardAxes,"XScale",hazardAxis.xScale, ...
    "XLim",hazardAxis.xLim,"XTick",hazardAxis.xTicks, ...
    "XTickLabel",hazardAxis.xTickLabels);

% Saving
if save_option
    output_directory = 'figures';

    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end

    savefig(fig, fullfile(output_directory, ...
        'canonical_parameters.fig'));

    exportgraphics(fig, fullfile(output_directory, ...
        'canonical_parameters.png'), ...
        'Resolution', 300);

    exportgraphics(fig, fullfile(output_directory, ...
        'canonical_parameters.pdf'), ...
        'ContentType', 'image','Resolution', 300);

    disp('Done');
end

%% Local functions
function sweep = makeSweep(field,values,symbol)
sweep.field = field;
sweep.values = values;
sweep.symbol = symbol;
end

function settings = makeTimeAxisSettings(scale, ...
    linearLim,linearTicks,linearTickLabels, ...
    logLim,logTicks,logTickLabels)
settings.xScale = scale;

if scale == "linear"
    settings.xLim = linearLim;
    settings.xTicks = linearTicks;
    settings.xTickLabels = linearTickLabels;
else
    settings.xLim = logLim;
    settings.xTicks = logTicks;
    settings.xTickLabels = logTickLabels;
end
end

function validateTimeAxisSettings(settings,timeEvalLimits,name)
if ~isscalar(settings.xScale) || ...
        ~any(settings.xScale == ["linear","log"])
    error(name+".xScale must be ""linear"" or ""log"".");
end
if numel(settings.xLim) ~= 2 || ...
        settings.xLim(2) <= settings.xLim(1) || ...
        settings.xLim(1) < timeEvalLimits(1) || ...
        settings.xLim(2) > timeEvalLimits(2)
    error(name+".xLim must be increasing and lie within timeEvalLimits.");
end
if any(settings.xTicks < settings.xLim(1)) || ...
        any(settings.xTicks > settings.xLim(2))
    error(name+".xTicks must lie inside its x-limits.");
end
if numel(settings.xTicks) ~= numel(settings.xTickLabels)
    error(name+" must have one label for every x-tick.");
end
if settings.xScale == "log" && ...
        (settings.xLim(1) <= 0 || any(settings.xTicks <= 0))
    error(name+" requires positive limits and ticks for logarithmic scale.");
end
end

function out = bifurcationData(p,z)
% For the unidirectional canonical cycle,
%   <T_FF> = 1/gamma0 + n/gamma + 1/A(z),
% and p_F = 1/(gamma0*<T_FF>).
[A,dAdz] = activationRate(z,p);
denominator = 1 + p.n*p.gamma0/p.gamma + p.gamma0./A;
pF = 1./denominator;

% Analytic derivative of p_F(z), followed by u(z)=z-w*p_F(z).
dpFdz = p.gamma0.*dAdz.*pF.^2./A.^2;
u = z - p.w*pF;
dudz = 1 - p.w*dpFdz;

% Saddle nodes occur at du/dz=0. Linear interpolation is adequate on the
% dense z grid used above.
idx = find(dudz(1:end-1).*dudz(2:end) < 0);
foldZ = zeros(size(idx));
foldU = zeros(size(idx));
foldP = zeros(size(idx));
for k = 1:numel(idx)
    i = idx(k);
    foldZ(k) = z(i) - dudz(i)*(z(i+1)-z(i))/(dudz(i+1)-dudz(i));
    foldU(k) = interp1(z,u,foldZ(k),'linear');
    foldP(k) = interp1(z,pF,foldZ(k),'linear');
end

out.z = z;
out.A = A;
out.pF = pF;
out.u = u;
out.dudz = dudz;
out.foldZ = foldZ;
out.foldU = foldU;
out.foldP = foldP;
out.recoveryCeiling = 1/(1 + p.n*p.gamma0/p.gamma);
end

function [A,dAdz] = activationRate(z,p)
% signedPower keeps z^nu real and monotone when z<0 and nu is noninteger.
% Replace x by max(z,0).^p.nu here if a rectified-power convention is wanted.
x = sign(z).*abs(z).^p.nu;
if p.nu == 1
    dxdz = ones(size(z));
else
    dxdz = p.nu*abs(z).^(p.nu-1);
end

q = p.thetaA*(x-p.muA);
q = max(min(q,700),-700);       % Avoid floating-point overflow.
sigma = 1./(1+exp(-q));
A = p.Amax*sigma;
A = max(A,realmin('double'));
dAdz = p.thetaA*A.*(1-A/p.Amax).*dxdz;
end

function eq = selectStableEquilibrium(p,uTarget,z,branch)
b = bifurcationData(p,z);
g = b.u-uTarget;
idx = find(g(1:end-1).*g(2:end) <= 0);

zRoots = zeros(size(idx));
for k = 1:numel(idx)
    i = idx(k);
    if g(i+1) == g(i)
        zRoots(k) = z(i);
    else
        zRoots(k) = z(i)-g(i)*(z(i+1)-z(i))/(g(i+1)-g(i));
    end
end
zRoots = unique(zRoots);
dudzRoots = interp1(z,b.dudz,zRoots,'linear');
zStable = zRoots(dudzRoots > 0);

if isempty(zStable)
    error(['No stable equilibrium was found at uEval=%.4g. ', ...
        'Increase the zGrid range or change uEval.'],uTarget);
end

switch lower(branch)
    case 'high'
        zEq = max(zStable);
    case 'low'
        zEq = min(zStable);
    otherwise
        error('equilibriumBranch must be ''high'' or ''low''.');
end

[AEq,~] = activationRate(zEq,p);
pFEq = 1/(1 + p.n*p.gamma0/p.gamma + p.gamma0/AEq);
eq.z = zEq;
eq.u = uTarget;
eq.A = AEq;
eq.pF = pFEq;
end

function [pdf,hazard] = returnTimeStatistics(p,AEq,t)
% Transient phases are [F,C1,...,Cn,Q]. Absorption is the next Q -> F jump.
rates = [p.gamma0,repmat(p.gamma,1,p.n),AEq];
nPhases = numel(rates);
S = diag(-rates);
for k = 1:nPhases-1
    S(k,k+1) = rates(k);
end
exitRates = -S*ones(nPhases,1);
alpha = [1,zeros(1,nPhases-1)];

pdf = zeros(size(t));
survival = zeros(size(t));
for k = 1:numel(t)
    transientProbability = alpha*expm(S*t(k));
    survival(k) = max(real(transientProbability*ones(nPhases,1)),0);
    pdf(k) = max(real(transientProbability*exitRates),0);
end

hazard = pdf./max(survival,eps);
end

function plotBifurcationPanel(ax,curve,colors,uLimits,uEval)
hold(ax,'on');
nCurves = numel(curve);

if curve(1).vary == "$\mathcal{A}_{\max}$"
    yMax = 0.12;
elseif curve(1).vary == "$\gamma$"
    yMax = 0.22;
else
    yMax = 0.1;
end

h = gobjects(1,nCurves);
for j = 1:nCurves
    b = curve(j).bif;
    stableP = b.pF;
    stableP(b.dudz <= 0) = NaN;
    unstableP = b.pF;
    unstableP(b.dudz > 0) = NaN;

    h(j) = plot(ax,b.u,stableP,'Color',colors(j,:), ...
        'LineWidth',1.5);
    % plot(ax,b.u,unstableP,':','Color',colors(j,:), ...
    %     'LineWidth',1.5,'HandleVisibility','off');

    % Filled circles identify the equilibria used for columns 2 and 3.
    xline(ax,uEval,'LineStyle','-','LineWidth',0.75,'Color',[0.8 0 0]);
    plot(ax,uEval,curve(j).eq.pF,'o','Color',colors(j,:), ...
        'MarkerFaceColor',colors(j,:),'MarkerSize',5, ...
        'HandleVisibility','off');
    xlim(ax,uLimits);
    ylim(ax,[0,yMax]);

    totalLine = plotVisualDots(ax,b.u,unstableP,colors(j,:),8,6);
end

xlabel(ax,'$u$','Interpreter','latex');
ylabel(ax,'$p_F^*$','Interpreter','latex');
formatAxis(ax);
leg = legend(ax,h,{curve.label},'Interpreter','latex','Location','west', ...
    'FontSize',10,'Box','on');
leg.Title.String = curve.vary;
leg.Title.FontSize = 11;
leg.Title.Interpreter = 'latex';
end

function formatAxis(ax)
grid(ax,'on');
box(ax,'on');
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 10;
ax.LineWidth = 0.5;
end

function addRowLabel(ax,rowNumber)
letter = char(double('a')+rowNumber-1);
text(ax,-0.125,1.2,sprintf('$(%c)$',letter), ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','top','Interpreter','latex', ...
    'FontSize',12,'FontWeight','bold');
end

function lineHandle = plotVisualDots(ax,x,y,color,spacingPixels,dotSize)

x = x(:);
y = y(:);

drawnow;

xLimits = xlim(ax);
yLimits = ylim(ax);

valid = isfinite(x) & isfinite(y) & ...
    x >= xLimits(1) & x <= xLimits(2) & ...
    y >= yLimits(1) & y <= yLimits(2);

if ax.XScale == "log"
    valid = valid & x > 0;
end

if ax.YScale == "log"
    valid = valid & y > 0;
end

x = x(valid);
y = y(valid);

if isempty(x)
    lineHandle = plot(ax,NaN,NaN, ...
        "LineStyle","none", ...
        "Marker",".", ...
        "MarkerSize",dotSize, ...
        "Color",color);
    return
end

xDisplay = axisForwardTransform(x,ax.XScale);
yDisplay = axisForwardTransform(y,ax.YScale);

xDisplay = xDisplay(:);
yDisplay = yDisplay(:);

xLimitDisplay = axisForwardTransform(xLimits,ax.XScale);
yLimitDisplay = axisForwardTransform(yLimits,ax.YScale);

axisPixels = getpixelposition(ax,true);

xPixels = axisPixels(3)*(xDisplay-xLimitDisplay(1))/diff(xLimitDisplay);
yPixels = axisPixels(4)*(yDisplay-yLimitDisplay(1))/diff(yLimitDisplay);

segmentLength = hypot(diff(xPixels),diff(yPixels));
arcLength = [0;cumsum(segmentLength(:))];

[arcLength,uniqueIndices] = unique(arcLength,"stable");

xDisplay = xDisplay(uniqueIndices);
yDisplay = yDisplay(uniqueIndices);

if isscalar(arcLength) || arcLength(end) == 0
    xDotsDisplay = xDisplay(1);
    yDotsDisplay = yDisplay(1);
else
    dotLocations = (0:spacingPixels:arcLength(end)).';

    if dotLocations(end) < arcLength(end)
        dotLocations(end+1,1) = arcLength(end);
    end

    xDotsDisplay = interp1( ...
        arcLength,xDisplay,dotLocations,"linear");

    yDotsDisplay = interp1( ...
        arcLength,yDisplay,dotLocations,"linear");
end

if ax.XScale == "log"
    xDots = 10.^xDotsDisplay;
else
    xDots = xDotsDisplay;
end

if ax.YScale == "log"
    yDots = 10.^yDotsDisplay;
else
    yDots = yDotsDisplay;
end

lineHandle = plot(ax,xDots,yDots, ...
    "LineStyle","none", ...
    "Marker",".", ...
    "MarkerSize",dotSize, ...
    "Color",color);

end

function transformed = axisForwardTransform(values,scale)

if scale == "log"
    transformed = log10(values);
else
    transformed = values;
end

end
