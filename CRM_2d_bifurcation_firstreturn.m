%% Three-pathway model: 2-D bifurcation map and return-statistic decomposition
% Canonical pathway:
%   F -> C1 -> ... -> Cn -> Q_C -> F
% Recurrence pathway:
%   C_k -> F
% Sequestration pathway (constant entry rate kappa_l):
%   C_l -> Y1 -> ... -> Ym -> Q_S -> F
%
% Q_C and Q_S are provenance-resolved copies of the same physical Q state.
% This augmentation is required to distinguish canonical and sequestration
% returns after the pathways merge. The physical Q occupancy is Q_C + Q_S.

clear; close all; clc;
save_option = 1;

%% User controls
% Canonical pathway
p.n       = 10;
p.gamma0  = 1;
p.gamma   = 1;
p.Amax    = 1;
p.muA     = 3.5;
p.thetaA  = 2.5;
p.w       = 8;
p.nu      = 1;

% Recurrence pathway C_k -> F, with rate beta_k P(z)
p.k       = 2;
p.Pmax    = 1;
p.muP     = 2.5;
p.thetaP  = 3.5;

% Sequestration pathway C_l -> Y1 -> ... -> Ym -> Q_S. The transition
% C_l -> Y1 occurs at the constant rate kappa_l, independent of z.
p.l       = 1;
p.m       = 30;
p.alpha   = 1*ones(1,p.m);  % Y_j forward rates, including Y_m -> Q_S

% Fixed external drive for the (beta_k,kappa_l) bifurcation map
uFixed = 1.75;

% Parameter plane
betaValues  = linspace(0,12,241);
kappaValues = linspace(0,0.2,101);
zGrid       = linspace(-8,18,3501);

% Seven operating points labelled (b)-(h): [beta_k,kappa_l]
operatingPoints = [
    2,0;   % (b): monostable low
    4,0;   % (c): bistable
    6,0;   % (d): monostable high
    8,0.15;   % (e): monostable low
    8,0.10;   % (f): bistable
    8,0.05; %(g): monostable high
    8,0% (h): monostable high
];

% Label offsets [Delta beta_k, Delta kappa_l] in panel (a)
pointLabelOffsets = [
   -0.07,0.005;
   -0.07,0.005;
   0.07,0.005;
    0.07,0.005;
    0.07,0.005;
    0.07,0.005;
    0.07,0.005
];

% Dashed parameter slices corresponding to the horizontal and vertical
% operating-point series. Set showGuideLines = false to suppress them.
showGuideLines = true;
betaGuide       = 8;
kappaGuide      = 0;

% Quantity displayed in panels (b)-(h): "hazard" or "first-return".
panelQuantity = "first-return";

% Time and return-statistic axes
timeEvalLimits = [0,100];
timePlotLimits = [0,100];
nTime          = 100001;
returnPanelYLimits = [];  % Empty: common automatic limits for every FP panel
highActivityThreshold = 0.12;

% Figure appearance
showBistableShading = true;
bistableShadeColor  = [0.82,0.82,0.82];
bistableShadeAlpha  = 0.5;
lowSNColor          = [1,0.5,0];
highSNColor         = [0.5,0,1];
pointColor          = [0.8,0,0];

canonicalColor      = [0.3,0.8,0.3];
recurrenceColor     = [1.00,0.30,0.70];
sequestrationColor  = [0.00,0.4,0.70];
totalColor          = [0,0,0];
meanReturnColor     = [0.8,0,0];
componentLineWidth  = 1.5;
totalLineWidth      = 2;
meanReturnLineWidth = 1.25;

%% Validation
validateParameters(p);
if numel(betaValues) < 3 || any(diff(betaValues) <= 0)
    error("betaValues must contain at least three strictly increasing values.");
end
if numel(kappaValues) < 3 || any(diff(kappaValues) <= 0)
    error("kappaValues must contain at least three strictly increasing values.");
end
if size(operatingPoints,1) ~= 7 || size(operatingPoints,2) ~= 2
    error("operatingPoints must be a 7-by-2 matrix for panels (b)-(h).");
end
if ~isequal(size(pointLabelOffsets),[7,2])
    error("pointLabelOffsets must be a 7-by-2 matrix.");
end
if any(operatingPoints(:,1) < betaValues(1) | ...
        operatingPoints(:,1) > betaValues(end)) || ...
        any(operatingPoints(:,2) < kappaValues(1) | ...
        operatingPoints(:,2) > kappaValues(end))
    error("Every operating point must lie inside the parameter plane.");
end
if timeEvalLimits(1) ~= 0 || timeEvalLimits(2) <= 0
    error("timeEvalLimits must have the form [0,tMax], with tMax > 0.");
end
if timePlotLimits(1) < timeEvalLimits(1) || ...
        timePlotLimits(2) > timeEvalLimits(2) || ...
        timePlotLimits(2) <= timePlotLimits(1)
    error("timePlotLimits must be increasing and lie within timeEvalLimits.");
end
if nTime < 2 || nTime ~= round(nTime)
    error("nTime must be an integer greater than one.");
end
if ~isscalar(panelQuantity) || ...
        ~any(panelQuantity == ["hazard","first-return"])
    error("panelQuantity must be ""hazard"" or ""first-return"".");
end
validateOptionalLimits(returnPanelYLimits,"returnPanelYLimits");

%% Two-dimensional bifurcation map
[Agrid,Pgrid] = regulatedFunctions(zGrid,p);
map = scanParameterPlane(p,zGrid,Agrid,Pgrid, ...
    betaValues,kappaValues,uFixed);

%% Fixed points and return statistics at the marked parameter values
t = linspace(timeEvalLimits(1),timeEvalLimits(2),nTime);
lowEquilibria = cell(1,7);
highEquilibria = cell(1,7);
lowStatistics = cell(1,7);
highStatistics = cell(1,7);

for j = 1:7
    beta = operatingPoints(j,1);
    kappa = operatingPoints(j,2);
    curve = equilibriumCurve(p,zGrid,Agrid,Pgrid,beta,kappa);
    stableZ = stableRootsAtDrive(zGrid,curve,uFixed);
    [lowEquilibria{j},highEquilibria{j}] = classifyEquilibria( ...
        p,zGrid,curve,stableZ,beta,kappa,highActivityThreshold);

    if ~isempty(lowEquilibria{j})
        lowStatistics{j} = returnStatistics(p,lowEquilibria{j},t);
    end
    if ~isempty(highEquilibria{j})
        highStatistics{j} = returnStatistics(p,highEquilibria{j},t);
    end
end

if isempty(returnPanelYLimits)
    returnPanelYLimits = automaticPanelLimits( ...
        [lowStatistics,highStatistics],panelQuantity);
end

%% Figure layout
fig = figure("Color","w","Units","pixels","Position",[20,120,1200,720]);
tl = tiledlayout(fig,41,51,"TileSpacing","tight","Padding","tight");

% Panel (a): 2-D bifurcation map spanning the upper-left region.
% Rows 1:20 and columns 1:30 are the fivefold-scaled version of the
% original 4-by-6 allocation. Row 21 is intentionally left empty to
% separate this map from panels (b)-(d), and column 31 is intentionally
% left empty to separate it from panels (e)-(h).
axA = nexttile(tl,1,[20,30]);
hold(axA,"on");

bistableMask = map.stableCount >= 2;
if showBistableShading
    shade = imagesc(axA,betaValues,kappaValues,double(bistableMask));
    shade.AlphaData = bistableShadeAlpha*double(bistableMask);
    colormap(axA,[1,1,1;bistableShadeColor]);
end
set(axA,"YDir","normal");

boundaryHandles = gobjects(0);
boundaryLabels = strings(0,1);

[contourLow,hLowSN] = contour(axA,betaValues,kappaValues, ...
    map.lowStateFoldU,[uFixed,uFixed],"-","Color",lowSNColor, ...
    "LineWidth",2);
if isempty(contourLow)
    delete(hLowSN);
else
    boundaryHandles(end+1) = hLowSN;
    boundaryLabels(end+1) = "Low-activity SN";
end

[contourHigh,hHighSN] = contour(axA,betaValues,kappaValues, ...
    map.highStateFoldU,[uFixed,uFixed],"-","Color",highSNColor, ...
    "LineWidth",2);
if isempty(contourHigh)
    delete(hHighSN);
else
    boundaryHandles(end+1) = hHighSN;
    boundaryLabels(end+1) = "High-activity SN";
end

if showGuideLines
    xline(axA,betaGuide,"--","Color",[0.35,0.35,0.35], ...
        "LineWidth",1.5,"HandleVisibility","off");
    yline(axA,kappaGuide,"--","Color",[0.35,0.35,0.35], ...
        "LineWidth",1.5,"HandleVisibility","off");
end

pointLetters = ["b","c","d","e","f","g","h"];
for j = 1:7
    beta = operatingPoints(j,1);
    kappa = operatingPoints(j,2);
    if pointLabelOffsets(j,1) < 0
        labelHorizontalAlignment = "right";
    else
        labelHorizontalAlignment = "left";
    end
    plot(axA,beta,kappa,"s","Color",pointColor, ...
        "MarkerFaceColor",pointColor,"MarkerSize",7,"LineWidth",1.0, ...
        "HandleVisibility","off");
    text(axA,beta+pointLabelOffsets(j,1), ...
        kappa+pointLabelOffsets(j,2),"$("+pointLetters(j)+")$", ...
        "Interpreter","latex","Color",pointColor,"FontSize",12, ...
        "HorizontalAlignment",labelHorizontalAlignment, ...
        "VerticalAlignment","bottom");
end

set(axA,'YLim',[kappaValues(1),kappaValues(end)],'YTick',0:0.05:0.15);
xlim(axA,[betaValues(1),betaValues(end)]);
xlabel(axA,"$\beta_k$","Interpreter","latex");
ylabel(axA,"$\kappa_l$","Interpreter","latex");
grid(axA,"on");
box(axA,"on");
axA.TickLabelInterpreter = "latex";
axA.LineWidth = 0.5;
axA.FontSize = 12;
text(axA,-0.05,0.985,"(a)", ...
    "Units","normalized","HorizontalAlignment","left", ...
    "VerticalAlignment","top","Interpreter","latex", ...
    "FontSize",13,"FontWeight","bold");

if ~isempty(boundaryHandles)
    legend(axA,boundaryHandles,boundaryLabels,"Interpreter","latex", ...
        "Location","NorthWest","FontSize",12,"Box","on");
end

% Return-statistic panel placement reproduces the reference structure: (b)-(d) are
% below panel (a), and (e)-(h) run down its right side. Each labelled panel
% consists of a low-FP axis and a high-FP axis. Every original tile is
% expanded to a 5-by-5 block. The lower-left panels begin at row 22, below
% the empty spacer row. The right-hand panels begin at column 32, beyond the
% empty spacer column, and end at row 40 so that row 41 is empty there.
lowTiles  = [1072,1082,1092, 32,542,1052,1562];
highTiles = [1582,1592,1602, 42,552,1062,1572];
lowAxes = gobjects(1,7);
highAxes = gobjects(1,7);
legendHandles = gobjects(0);

for j = 1:7
    lowAxes(j) = nexttile(tl,lowTiles(j),[10,10]);
    handles = plotReturnPanel(lowAxes(j),t,lowStatistics{j}, ...
        lowEquilibria{j},"low", ...
        canonicalColor,recurrenceColor,sequestrationColor,totalColor, ...
        meanReturnColor,componentLineWidth,totalLineWidth, ...
        meanReturnLineWidth,timePlotLimits,returnPanelYLimits,panelQuantity);
    if j <= 3
        addOuterPanelLabel(lowAxes(j),pointLetters(j),"above");
    else
        addOuterPanelLabel(lowAxes(j),pointLetters(j),"left");
    end
    if isempty(legendHandles) && ~isempty(handles)
        legendHandles = handles;
    end

    highAxes(j) = nexttile(tl,highTiles(j),[10,10]);
    handles = plotReturnPanel(highAxes(j),t,highStatistics{j}, ...
        highEquilibria{j},"high", ...
        canonicalColor,recurrenceColor,sequestrationColor,totalColor, ...
        meanReturnColor,componentLineWidth,totalLineWidth, ...
        meanReturnLineWidth,timePlotLimits,returnPanelYLimits,panelQuantity);
    if isempty(legendHandles) && ~isempty(handles)
        legendHandles = handles;
    end
end

allReturnAxes = [lowAxes,highAxes];
linkaxes(allReturnAxes,"xy");


if ~isempty(legendHandles)
    if panelQuantity == "hazard"
        legendLabels = ["$w^C(t) h^C(t)$","$w^R(t) h^R(t)$","$w^S(t) h^S(t)$","$h(t)$"];
    else
        legendLabels = ["$f^C(t)$","$f^R(t)$","$f^S(t)$","$f(t)$"];
    end

    % Use the empty high-activity axis in panel (e) for the shared pathway
    % legend. Keep the axes and its state annotation, but move the
    % "not present" label upward to leave room for a vertical legend.
    legendAx = highAxes(4);
    notPresentLabel = findobj(legendAx,"Tag","NotPresentLabel");
    if ~isempty(notPresentLabel)
        set(notPresentLabel,"Position",[0.5,0.75,0]);
    end
    hold(legendAx,"on");
    legendDisplayHandles = [ ...
        plot(legendAx,nan,nan,"-","Color",canonicalColor, ...
            "LineWidth",componentLineWidth), ...
        plot(legendAx,nan,nan,"-","Color",recurrenceColor, ...
            "LineWidth",componentLineWidth), ...
        plot(legendAx,nan,nan,"-","Color",sequestrationColor, ...
            "LineWidth",componentLineWidth), ...
        plot(legendAx,nan,nan,":","Color",totalColor, ...
            "LineWidth",totalLineWidth)];
    legend(legendAx,legendDisplayHandles, ...
        legendLabels, ...
        "Interpreter","latex","Orientation","vertical", ...
        "NumColumns",1,"FontSize",11,"Box","on", ...
        "Location",[0.87 0.825 0.1 0.11]);
end

% Saving
if save_option
    output_directory = 'figures';

    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end


    savefig(fig, fullfile(output_directory, ...
        "treatment_2Dbif_"+panelQuantity+".fig"));

    exportgraphics(fig, fullfile(output_directory, ...
        "treatment_2Dbif_"+panelQuantity+".png"), ...
        'Resolution', 300);

    exportgraphics(fig, fullfile(output_directory, ...
        "treatment_2Dbif_"+panelQuantity+".pdf"), ...
        'ContentType', 'image','Resolution', 300);

    disp('Done');
end

%% Local functions
function map = scanParameterPlane(p,z,A,P,betaValues,kappaValues,uFixed)
nBeta = numel(betaValues);
nKappa = numel(kappaValues);
stableCount = zeros(nKappa,nBeta);
lowStateFoldU = nan(nKappa,nBeta);
highStateFoldU = nan(nKappa,nBeta);

for iy = 1:nKappa
    kappa = kappaValues(iy);
    for ix = 1:nBeta
        beta = betaValues(ix);
        curve = equilibriumCurve(p,z,A,P,beta,kappa);
        stableCount(iy,ix) = numel(stableRootsAtDrive(z,curve,uFixed));

        idx = find(curve.dudz(1:end-1).*curve.dudz(2:end) < 0);
        if numel(idx) >= 2
            foldZ = zeros(size(idx));
            for j = 1:numel(idx)
                i = idx(j);
                foldZ(j) = z(i)-curve.dudz(i)*(z(i+1)-z(i))/ ...
                    (curve.dudz(i+1)-curve.dudz(i));
            end
            foldZ = sort(foldZ);
            lowStateFoldU(iy,ix) = interp1(z,curve.u,foldZ(1),"linear");
            highStateFoldU(iy,ix) = interp1(z,curve.u,foldZ(end),"linear");
        end
    end
end

map.stableCount = stableCount;
map.lowStateFoldU = lowStateFoldU;
map.highStateFoldU = highStateFoldU;
end

function curve = equilibriumCurve(p,z,A,P,beta,kappa)
recurrenceRate = beta.*P;
sequestrationRate = kappa*ones(size(z));
arrival = ones(size(z));
meanTransientTime = zeros(size(z));
rhoR = zeros(size(z));
rhoS = zeros(size(z));

for i = 1:p.n
    totalRate = p.gamma*ones(size(z));
    if i == p.k
        totalRate = totalRate+recurrenceRate;
    end
    if i == p.l
        totalRate = totalRate+sequestrationRate;
    end

    meanTransientTime = meanTransientTime+arrival./totalRate;
    if i == p.k
        rhoR = rhoR+arrival.*recurrenceRate./totalRate;
    end
    if i == p.l
        rhoS = rhoS+arrival.*sequestrationRate./totalRate;
    end
    arrival = arrival.*p.gamma./totalRate;
end

rhoC = arrival;
meanTransientTime = meanTransientTime ...
    + rhoS*sum(1./p.alpha) ...
    + (rhoC+rhoS)./A;
meanCycleTime = 1/p.gamma0+meanTransientTime;
pF = (1/p.gamma0)./meanCycleTime;

curve.pF = pF;
curve.u = z-p.w*pF;
curve.dudz = gradient(curve.u,z);
curve.A = A;
curve.recurrenceRate = recurrenceRate;
curve.sequestrationRate = sequestrationRate;
curve.rhoC = rhoC;
curve.rhoR = rhoR;
curve.rhoS = rhoS;
end

function stableZ = stableRootsAtDrive(z,curve,uTarget)
g = curve.u-uTarget;
idx = find(g(1:end-1).*g(2:end) <= 0);
roots = zeros(size(idx));
for j = 1:numel(idx)
    i = idx(j);
    denominator = g(i+1)-g(i);
    if denominator == 0
        roots(j) = z(i);
    else
        roots(j) = z(i)-g(i)*(z(i+1)-z(i))/denominator;
    end
end
roots = unique(roots);
if isempty(roots)
    stableZ = [];
    return;
end
dudzRoots = interp1(z,curve.dudz,roots,"linear");
stableZ = sort(roots(dudzRoots > 0));
end

function [lowEq,highEq] = classifyEquilibria( ...
    p,z,curve,stableZ,beta,kappa,highActivityThreshold)
lowEq = [];
highEq = [];
if isempty(stableZ)
    return;
end

if numel(stableZ) >= 2
    lowEq = makeEquilibrium(p,z,curve,stableZ(1),beta,kappa);
    highEq = makeEquilibrium(p,z,curve,stableZ(end),beta,kappa);
else
    eq = makeEquilibrium(p,z,curve,stableZ(1),beta,kappa);
    if eq.pF >= highActivityThreshold
        highEq = eq;
    else
        lowEq = eq;
    end
end
end

function eq = makeEquilibrium(p,z,curve,zEq,beta,kappa)
eq.beta = beta;
eq.kappa = kappa;
eq.z = zEq;
eq.pF = interp1(z,curve.pF,zEq,"linear");
eq.A = interp1(z,curve.A,zEq,"linear");
eq.recurrenceRate = interp1(z,curve.recurrenceRate,zEq,"linear");
eq.sequestrationRate = interp1(z,curve.sequestrationRate,zEq,"linear");
eq.rhoC = interp1(z,curve.rhoC,zEq,"linear");
eq.rhoR = interp1(z,curve.rhoR,zEq,"linear");
eq.rhoS = interp1(z,curve.rhoS,zEq,"linear");

rhoSum = eq.rhoC+eq.rhoR+eq.rhoS;
eq.rhoC = eq.rhoC/rhoSum;
eq.rhoR = eq.rhoR/rhoSum;
eq.rhoS = eq.rhoS/rhoSum;
end

function out = returnStatistics(p,eq,t)
% Transient states: [F,C_1,...,C_n,Y_1,...,Y_m,Q_C,Q_S].
nStates = p.n+p.m+3;
iF = 1;
iC = 2:(p.n+1);
iY = (p.n+2):(p.n+p.m+1);
iQC = p.n+p.m+2;
iQS = p.n+p.m+3;

S = zeros(nStates);
S(iF,iF) = -p.gamma0;
S(iF,iC(1)) = p.gamma0;

exitCanonical = zeros(nStates,1);
exitRecurrence = zeros(nStates,1);
exitSequestration = zeros(nStates,1);

for i = 1:p.n
    state = iC(i);
    totalRate = p.gamma;
    if i == p.k
        totalRate = totalRate+eq.recurrenceRate;
        exitRecurrence(state) = eq.recurrenceRate;
    end
    if i == p.l
        totalRate = totalRate+eq.sequestrationRate;
        S(state,iY(1)) = eq.sequestrationRate;
    end
    S(state,state) = -totalRate;
    if i < p.n
        S(state,iC(i+1)) = p.gamma;
    else
        S(state,iQC) = p.gamma;
    end
end

for j = 1:p.m
    S(iY(j),iY(j)) = -p.alpha(j);
    if j < p.m
        S(iY(j),iY(j+1)) = p.alpha(j);
    else
        S(iY(j),iQS) = p.alpha(j);
    end
end

S(iQC,iQC) = -eq.A;
S(iQS,iQS) = -eq.A;
exitCanonical(iQC) = eq.A;
exitSequestration(iQS) = eq.A;

probability = zeros(1,nStates);
probability(iF) = 1;
meanReturnTime = -probability*(S\ones(nStates,1));
dt = t(2)-t(1);
transition = expm(S*dt);

hCanonical = zeros(size(t));
hRecurrence = zeros(size(t));
hSequestration = zeros(size(t));
fCanonical = zeros(size(t));
fRecurrence = zeros(size(t));
fSequestration = zeros(size(t));
survival = zeros(size(t));

for it = 1:numel(t)
    survival(it) = max(real(sum(probability)),0);
    fCanonical(it) = max(real(probability*exitCanonical),0);
    fRecurrence(it) = max(real(probability*exitRecurrence),0);
    fSequestration(it) = max(real(probability*exitSequestration),0);
    if survival(it) > realmin("double")
        hCanonical(it) = fCanonical(it)/survival(it);
        hRecurrence(it) = fRecurrence(it)/survival(it);
        hSequestration(it) = fSequestration(it)/survival(it);
    end
    probability = probability*transition;
end

out.canonical = hCanonical;
out.recurrence = hRecurrence;
out.sequestration = hSequestration;
out.total = hCanonical+hRecurrence+hSequestration;
out.survival = survival;
out.fCanonical = fCanonical;
out.fRecurrence = fRecurrence;
out.fSequestration = fSequestration;
out.fTotal = out.fCanonical+out.fRecurrence+out.fSequestration;
out.meanReturnTime = meanReturnTime;
end

function handles = plotReturnPanel(ax,t,statistics,eq,branch, ...
    canonicalColor,recurrenceColor,sequestrationColor,totalColor, ...
    meanReturnColor,componentLineWidth,totalLineWidth, ...
    meanReturnLineWidth,xLimits,yLimits,panelQuantity)
hold(ax,"on");

if branch == "low"
    activityLabel = "Low-activity FP";
else
    activityLabel = "High-activity FP";
end
leftAnnotation = activityLabel;
rhoAnnotation = strings(0,1);

if isempty(eq)
    handles = gobjects(0);
    text(ax,0.5,0.53,"Not present", ...
        "Units","normalized","HorizontalAlignment","center", ...
        "VerticalAlignment","middle","Interpreter","latex","FontSize",11, ...
        "Tag","NotPresentLabel");
else
    if panelQuantity == "hazard"
        canonical = statistics.canonical;
        recurrence = statistics.recurrence;
        sequestration = statistics.sequestration;
        total = statistics.total;
    else
        canonical = statistics.fCanonical;
        recurrence = statistics.fRecurrence;
        sequestration = statistics.fSequestration;
        total = statistics.fTotal;
    end

    hCan = plot(ax,t,canonical,"-","Color",canonicalColor, ...
        "LineWidth",componentLineWidth);
    hRec = plot(ax,t,recurrence,"-","Color",recurrenceColor, ...
        "LineWidth",componentLineWidth);
    hSeq = plot(ax,t,sequestration,"-","Color",sequestrationColor, ...
        "LineWidth",componentLineWidth);
    hTotal = plot(ax,t,total,":","Color",totalColor, ...
        "LineWidth",totalLineWidth);
    scatter(ax,statistics.meanReturnTime,0.001,80,'Marker','*',"MarkerEdgeColor",meanReturnColor, ...
        'LineWidth',1.5,"HandleVisibility","off");
    handles = [hCan,hRec,hSeq,hTotal];

    leftAnnotation = [
        activityLabel
        "$z^* = "+compose("%.3f",eq.z)+",\ p_F^* = "+compose("%.3f",eq.pF)+"$"
    ];
    rhoAnnotation = [
        "$\rho^C = "+compose("%.3f",eq.rhoC)+"$"
        "$\rho^R = "+compose("%.3f",eq.rhoR)+"$"
        "$\rho^S = "+compose("%.3f",eq.rhoS)+"$"
    ];
end

text(ax,0.08,0.96,leftAnnotation,"Units","normalized", ...
    "HorizontalAlignment","left","VerticalAlignment","top", ...
    "Interpreter","latex","FontSize",10.5);
if ~isempty(rhoAnnotation)
    text(ax,0.94,0.58,rhoAnnotation,"Units","normalized", ...
        "HorizontalAlignment","right","VerticalAlignment","top", ...
        "Interpreter","latex","FontSize",9.5);
end
% set(ax,'XLim',[0.1 200],'XScale','log','XTick',[0.1 1 10 100]);
set(ax,'XLim',[0 90],'XScale','linear','XTick',[0 30 60 90]);

set(ax,'YLim',[0.001 6],'YScale','log','YTick',[0.001 0.01 0.1 1]);
% set(ax,'YLim',[0 0.55],'YScale','linear','YTick',[0 0.25 0.5]);

xlabel(ax,"$t$","Interpreter","latex");
if panelQuantity == "hazard"
    ylabel(ax,"$h(t)$","Interpreter","latex");
else
    ylabel(ax,"$f(t)$","Interpreter","latex");
end
grid(ax,"on");
box(ax,"on");
ax.TickLabelInterpreter = "latex";
ax.FontSize = 10;
ax.LineWidth = 0.5;
end

function limits = automaticPanelLimits(statisticCells,panelQuantity)
maximum = 0;
for j = 1:numel(statisticCells)
    if ~isempty(statisticCells{j})
        if panelQuantity == "hazard"
            total = statisticCells{j}.total;
        else
            total = statisticCells{j}.fTotal;
        end
        maximum = max(maximum,max(total));
    end
end
limits = [0,max(1e-6,1.06*maximum)];
end

function [A,P] = regulatedFunctions(z,p)
xA = sign(z).*abs(z).^p.nu;
qA = max(min(p.thetaA*(xA-p.muA),700),-700);
A = p.Amax./(1+exp(-qA));
A = max(A,realmin("double"));

qP = max(min(p.thetaP*(z-p.muP),700),-700);
P = p.Pmax./(1+exp(-qP));
end

function validateParameters(p)
if p.n < 1 || p.n ~= round(p.n)
    error("p.n must be a positive integer.");
end
if p.k < 1 || p.k > p.n || p.k ~= round(p.k)
    error("p.k must be an integer between 1 and p.n.");
end
if p.l < 1 || p.l > p.n || p.l ~= round(p.l)
    error("p.l must be an integer between 1 and p.n.");
end
if p.m < 1 || p.m ~= round(p.m) || numel(p.alpha) ~= p.m
    error("p.m must be positive and numel(p.alpha) must equal p.m.");
end
if any(~isfinite(p.alpha)) || any(p.alpha <= 0)
    error("Every sequestration recovery rate in p.alpha must be positive.");
end

positiveFields = ["gamma0","gamma","Amax","thetaA","w","nu", ...
    "Pmax","thetaP"];
for j = 1:numel(positiveFields)
    value = p.(positiveFields(j));
    if ~isscalar(value) || ~isfinite(value) || value <= 0
        error("p."+positiveFields(j)+" must be a positive finite scalar.");
    end
end
end

function validateOptionalLimits(limits,name)
if ~isempty(limits) && ...
        (numel(limits) ~= 2 || limits(1) < 0 || limits(2) <= limits(1))
    error(name+" must be empty or a nonnegative increasing pair.");
end
end





function addOuterPanelLabel(ax,letter,position)
if position == "above"
    x = -0.17;
    y = 1.01;
    horizontalAlignment = "left";
    verticalAlignment = "bottom";
else
    x = -0.15;
    y = 1;
    horizontalAlignment = "right";
    verticalAlignment = "top";
end
text(ax,x,y,"("+letter+")", ...
    "Units","normalized","HorizontalAlignment",horizontalAlignment, ...
    "VerticalAlignment",verticalAlignment,"Interpreter","latex", ...
    "FontSize",12,"Clipping","off");
end
