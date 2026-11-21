%% gillespie_QFCY_feedback.m
% Exact Gillespie simulation of N coupled cells on the Q-F-C-Y CTMC with
% QUENCHED parameter heterogeneity and firing-occupancy feedback.
%
%   Canonical pathway      Q --A(z)--> F --g0--> C1 --g1--> ... --> Cn --gn--> Q
%   Recurrence pathway     C_k --beta_k * P_k(z)-->  F
%   Sequestration pathway  C_l --kappa_l * D_l(z)--> Y1 --a1--> ... --> Ym --am--> Q
%
%   Sigmoids:  S(z) = 1 / (1 + exp(-theta*(z - mu)))       (S_max = 1)
%
%   Input:     z(t) = w * p_F(t) + u(t),   p_F = N_F / N  (firing occupancy)
%
% Because z depends on the whole population, cells are simulated jointly:
% every event is one cell changing state, and z changes whenever a cell
% enters or leaves F.
%
% ALGORITHM: population-level SSA with thinning (exact).
%   Every sigmoid is <= 1, so the rate of each exit is bounded by its
%   prefactor NO MATTER WHAT z IS. The total bound B = sum over cells of the
%   prefactors of their current exits therefore only changes when a cell
%   jumps; nothing has to be recomputed when z changes (via p_F or u(t)).
%     1. draw a candidate time at rate B
%     2. pick a cell with probability proportional to its bound (Fenwick
%        tree, O(log N) per draw/update)
%     3. evaluate z = w*p_F + u(t) at the candidate time and accept an exit
%        with probability (true rate)/(bound), otherwise reject.
%   This is exact for any w (either sign) and any u(t), even discontinuous.
%
% QUENCHED HETEROGENEITY: each cell draws its parameters once at t = 0 and
% keeps them for the whole run (log-normal rates/steepnesses with given CV,
% Gaussian sigmoid midpoints with given SD). Set all CVs/SDs to 0 for the
% homogeneous model.
%
% MEAN-FIELD CHECK: each cell's state probability p_i(t) obeys its own master
% equation dp_i/dt = p_i G_i(z) for the SAME sampled parameters. Two versions:
%   closed loop  z = w*mean_i p_i,F + u(t)   -- the N -> infinity limit. The
%                model is bistable, so near the separatrix finite-N noise can
%                send the SSA to the other attractor; the closed loop cannot
%                follow that (it is deterministic).
%   driven       z = z_SSA(t) = w*p_F,SSA(t) + u(t), the input the cells in the
%                SSA actually received. Given the input path the cells are
%                independent, so this is the expected occupancy conditioned on
%                the realised feedback: it stays in whichever attractor the SSA
%                selected, and any remaining mismatch is single-cell noise.
%
% Requires MATLAB R2019b or newer (local functions in scripts, tiledlayout).

clear; clc; rng(205);
save_option = 1;

%% ======================= Model structure =================================
par.n = 10;          % canonical stages C1..Cn
par.m = 50;          % sequestration stages Y1..Ym
par.k = 2;          % recurrence exit stage(s)    C_k -> F   (scalar or vector)
par.l = 1;          % sequestration exit stage(s) C_l -> Y1  (scalar or vector)

%% ================= Population-mean parameter values ======================
par.gamma = [3, 3*ones(1, par.n)];     % [gamma_0 ... gamma_n]
par.alpha = 3 * ones(1, par.m);   % [alpha_1 ... alpha_m]
par.beta  = 6;                    % beta_k  (one per element of par.k)
par.kappa = 1;                    % kappa_l (one per element of par.l)

par.thA = 2;   par.muA = 3;      % activation     A(z)
par.thP = 3.5;   par.muP = 2.5;      % recurrence     P_k(z)
par.thD = 5.5;   par.muD = 8.0;      % sequestration  D_l(z)

%% ===================== Quenched heterogeneity ============================
% CVs (log-normal) for rates and steepnesses, SDs (Gaussian) for midpoints
het.gamma = 0.20;   het.alpha = 0.20;   het.beta = 0.30;   het.kappa = 0.30;
het.thA   = 0.10;   het.thP   = 0.10;   het.thD  = 0.10;
het.muA   = 0.05;   het.muP   = 0.05;   het.muD  = 0.15;

%% ===================== Feedback and external input =======================
sim.w = 8;                                      % coupling weight
sim.u = @(t) 1.5 + 10 * (t >= 50 & t <= 51) + 5 * (t > 51); % external input u(t) (or a constant)

%% ======================= Simulation settings =============================
sim.N         = 2000;   % number of cells
sim.T         = 150;    % final time
sim.dt        = 0.25;   % output grid spacing
sim.s0        = ones(sim.N,1);      % initial state (1 = Q, 2 = F, 3.. = C1.., then Y1..);
                        % scalar or N-vector
sim.s0 = round(10*rand(sim.N,1)) + sim.s0;
sim.nTrace    = 4;      % cells whose full event history is stored
sim.meanField = 'both'; % 'driven' (by SSA feedback), 'closed' (N -> inf), 'both', or 'none'
sim.meanFieldHom = true; % also the homogeneous mean field: one cell with the population-mean
                        % parameters (62 equations), same initial state distribution

%% ===================== Spike-train analysis settings =====================
% A "spike" is any entry into F: activation Q->F or recurrence C_k->F.
ana.win     = 3;        % counting window (the 't' interval)
ana.step    = 0.3;      % window stride; set ana.step = ana.win for
                        % non-overlapping windows (ana.win/ana.step must be an integer)
ana.rasterN = 400;      % neurons shown in the raster (random subset); Inf = all
ana.rasterFrac   = [0.5 0.3 0.2];   % target share of raster rows per pathway class:
                                    % [recurrence R, canonical C, sequestered S]
ana.rasterRatio  = 3;               % raster display only: class R needs >= this many
                                    % recurrences per canonical return
ana.alignN    = 10;     % raster bottom: groups of this many neurons whose silent gaps start
ana.alignStep = 10;     % at tPulse, tPulse + alignStep, ... (bottom group = right after the pulse)
ana.alignTol  = 2;      % a gap belongs to anchor a if its last spike falls in [a, a + alignTol)
ana.alignJitter  = 4;   % anchors after the first are shifted by a random amount in +/- this
ana.alignShuffle = true;% stack the later groups in random order (bottom group stays the pulse)
ana.rasterBlocks = true;            % true: rows grouped R | C | S, each block sorted by
                                    % spike count; false: all rows sorted by spike count
ana.minISI  = 20;       % min. number of ISIs in a window to report the ISI CV
ana.rateTol = 0.33;      % x in the rate-change panel: > +x increasing, < -x decreasing,
                        % within +/-x same (0.2 = 20%)
ana.rateRef = 50 + [-40 0];  % baseline interval: each neuron's own mean rate over the
                        % 20 time units before the pulse (t = 50), i.e. [30, 50)
ana.rateSkipSilent = false; % don't count neurons with zero rate in both (no measurable change)
ana.tPulse  = 50;       % pulse onset: rate-change fractions are only shown for t > tPulse

%% ============================== Run ======================================
tic;
out = simulateCoupledQFCY(par, het, sim);
fprintf('SSA: %d cells, T = %g, %d events (acceptance %.1f%%), %.2f s\n', ...
        sim.N, sim.T, out.nAccepted, 100*out.nAccepted/out.nCandidates, toc);
for q = 1:numel(out.typeNames)
    fprintf('  %-30s %8.2f per cell\n', out.typeNames{q}, mean(out.counts(:, q)));
end

%% ============================== Plots ====================================
% All text (titles, labels, legends, tick labels) uses the LaTeX interpreter.
% Set as figure-level defaults so the MATLAB session defaults are untouched.
latexDefaults = {'DefaultTextInterpreter',          'latex', ...
                 'DefaultAxesTickLabelInterpreter', 'latex', ...
                 'DefaultLegendInterpreter',        'latex'};

t  = out.t;
iC = 2 + (1:par.n);   iY = 2 + par.n + (1:par.m);
group = @(Fr) [Fr(1,:); Fr(2,:); sum(Fr(iC,:), 1); sum(Fr(iY,:), 1)];
cols  = [0 0.55 0; 0.85 0.3 0; 0 0 0; 0 0.45 0.95];
names = {'$Q$', '$F$', '$C_1 \ldots C_n$', '$Y_1 \ldots Y_m$'};

figure('Color', 'w', 'Position', [60 60 1200 780], latexDefaults{:});
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- occupancy --------------------------------------------------------------
nexttile; hold on; box on;
G = group(out.frac);
for q = 1:4
    plot(t, G(q,:), '-', 'Color', cols(q,:), 'LineWidth', 1.5, 'DisplayName', names{q});
end
mfSets = {};                                   % {field, line style, colour, width, legend}
if isfield(out, 'fracMF')
    mfSets(end+1, :) = {'fracMF',  '-.', [0.6 0.6 0.6], 1.0, 'closed-loop mean field ($N \rightarrow \infty$)'};
end
if isfield(out, 'fracDMF')
    mfSets(end+1, :) = {'fracDMF', ':',  [0 0 0],       1.3, 'mean field driven by SSA $z(t)$'};
end
for k = 1:size(mfSets, 1)
    GM = group(out.(mfSets{k, 1}));
    for q = 1:4
        h = plot(t, GM(q,:), mfSets{k, 2}, 'Color', mfSets{k, 3}, 'LineWidth', mfSets{k, 4});
        if q == 1, h.DisplayName = mfSets{k, 5};
        else,      h.HandleVisibility = 'off'; end
    end
end
xlabel('time'); ylabel('fraction of cells'); legend('Location', 'northwest', 'Orientation', 'vertical');
title('Population occupancy');  

% --- input z(t) --------------------------------------------------------------
nexttile; hold on; box on;
plot(t, out.z, '-', 'Color', [0.85 0.3 0], 'LineWidth', 1.4, 'DisplayName', '$z = w\,p_F + u$ (SSA)');
if isfield(out, 'zDMF')
    plot(t, out.zDMF, 'k:', 'LineWidth', 1.3, 'DisplayName', '$z$ driven mean field');
end
if isfield(out, 'zMF')
    plot(t, out.zMF, '-.', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0, ...
         'DisplayName', '$z$ closed-loop mean field');
end
plot(t, out.u, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, 'DisplayName', '$u(t)$');
xlabel('time'); ylabel('$z$'); legend('Location', 'northwest', 'Orientation', 'vertical');
title(sprintf('Input with firing-occupancy feedback ($w = %g$)', sim.w));  

% --- individual cells ----------------------------------------------------------
nexttile; hold on; box on;
for i = 1:numel(out.traces)
    stairs(out.traces(i).t, out.traces(i).s + 0.1*(i-1), 'LineWidth', 1.1);
end
tickLab = strcat('$', regexprep(out.labels, '(\d+)', '_{$1}'), '$');   % C12 -> $C_{12}$
set(gca, 'YTick', 1:numel(out.labels), 'YTickLabel', tickLab, 'YDir', 'reverse');
ylim([0.5, numel(out.labels) + 0.5]); xlabel('time');
title('Individual cells (frozen parameters)'); 

% --- per-cell pathway usage -----------------------------------------------------
nexttile; hold on; box on;
histogram(out.counts(:, 4), 'BinMethod', 'integers', 'Normalization', 'probability', ...
          'FaceColor', [0.85 0.1 0.45], 'FaceAlpha', 0.6);
histogram(out.counts(:, 5), 'BinMethod', 'integers', 'Normalization', 'probability', ...
          'FaceColor', [0 0.45 0.95], 'FaceAlpha', 0.6);
legend({'recurrence $C_k \rightarrow F$', 'sequestration $C_l \rightarrow Y_1$'}, ...
       'Location', 'northwest', 'Orientation', 'vertical');
xlabel('events per cell in $[0, T]$'); ylabel('fraction of cells');
title('Cell-to-cell variability in pathway usage');  

%% ====================== Figures to save ===========================
%%% Figure 1
fig1 = figure('Color', 'w', 'Position', [60 60 600 350]);
% --- occupancy --------------------------------------------------------------
hold on; box on;
G = group(out.frac);
for q = 1:4
    plot(t, G(q,:), '-', 'Color', cols(q,:), 'LineWidth', 1.5, 'DisplayName', names{q});
end

% for k = 1:size(mfSets, 1)
%     GM = group(out.(mfSets{k, 1}));
%     for q = 1:4
%         h = plot(t, GM(q,:), ':', 'Color', cols(q,:), 'LineWidth', 2,'HandleVisibility','off');
%     end
% end
if isfield(out, 'fracHMF')                           % homogeneous mean field (mean parameters)
    GH = group(out.fracHMF);
    for q = 1:4
        plot(t, GH(q,:), ':', 'Color', cols(q,:), 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
% % line-style key (black, so it applies to every colour)
% if size(mfSets, 1) > 0
%     plot(NaN, NaN, 'k:', 'LineWidth', 2, 'DisplayName', 'Mean field (sampled parameters)');
% end
% if isfield(out, 'fracHMF')
%     plot(NaN, NaN, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Homogeneous mean field');
% end
xlabel("Time $t$","FontSize",12,"Interpreter","latex"); 
ylabel("$p_s^*$","FontSize",12,"Interpreter","latex");
legend('Location', 'best', 'Orientation', 'vertical',"FontSize",12,"Interpreter","latex");
set(gca,'XGrid','On','YGrid','On',"FontSize",12,'TickLabelInterpreter','latex');
title('Gillespie simulation vs Homogeneous mean-field limit',"FontSize",13,"Interpreter","latex");  

% Saving
if save_option
    output_directory = 'figures';

    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end
    

    savefig(fig1, fullfile(output_directory, ...
        "Gillespie_three_pathway.fig"));

    exportgraphics(fig1, fullfile(output_directory, ...
        "Gillespie_three_pathway.png"), ...
        'Resolution', 300);

    exportgraphics(fig1, fullfile(output_directory, ...
        "Gillespie_three_pathway.pdf"), ...
        'ContentType', 'image','Resolution', 300);

    disp('Done');
end


%%% Figure 2
S = spikeStatistics(out.spikes, sim.N, sim.T, ana.win, ana.step, ana.minISI);
fprintf('%d spikes recorded (%.1f per neuron)\n', numel(out.spikes.t), numel(out.spikes.t)/sim.N);
cls  = pathwayClass(out.counts, ana.rasterRatio);        % raster selection only: 1 R, 2 C, 3 S, 0 none
clsNames = {'R', 'C', 'S', 'other'};

fig2 = figure('Color', 'w', 'Position', [100 40 900 1100], latexDefaults{:});
tiledlayout(7, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
ax = gobjects(6, 1);

% --- raster -------------------------------------------------------------------
% Stratified random subset: each neuron is assigned a pathway class from its
% event counts (see pathwayClass), and rows are filled in the proportions
% ana.rasterFrac so that rare classes (sequestered cells) are always shown.
ax(1) = nexttile([2 1]); hold on; box on;
nShow = min(ana.rasterN, sim.N);
tot   = accumarray(out.spikes.cell, 1, [sim.N 1]);       % spikes per neuron
[sel, grp] = selectRasterCells(cls, tot, nShow, ana.rasterFrac, ana.rasterBlocks);
nShow = numel(sel);
fprintf('Raster rows (population):');
for g = 1:4
    fprintf('  %s %d (%d)', clsNames{g}, sum(grp == g), sum(cls == mod(g, 4)));
end
fprintf('\n');
% bottom of the raster: groups of neurons whose longest silence starts at the
% same time, one group per anchor (bottom group = right after the pulse, the next
% one alignStep later, ...). They replace rows of the random selection.
anchors = ana.tPulse : ana.alignStep : sim.T - ana.alignStep;
nA      = numel(anchors);
anchors(2:end) = anchors(2:end) + ana.alignJitter * (2 * rand(1, nA - 1) - 1);   % irregular times
anchors(2:end) = min(max(anchors(2:end), ana.tPulse + ana.alignTol), sim.T - ana.alignTol);
[alignIds, alignLen] = alignedSilences(out.spikes, sim.N, sim.T, anchors, ana.alignTol, ana.alignN);
stackOrder = 1:nA;                                        % bottom to top
if ana.alignShuffle, stackOrder = [1, 1 + randperm(nA - 1)]; end
nAlign  = sum(cellfun(@numel, alignIds));
selMain = sel(~ismember(sel, vertcat(alignIds{:})));
selMain = selMain(1:min(numel(selMain), nShow - nAlign));
sel     = [selMain(:); vertcat(alignIds{fliplr(stackOrder)})];   % pulse group at the very bottom
nShow   = numel(sel);
for q = 1:numel(anchors)
    fprintf('  aligned group t = %.1f: %d neurons, silent for %.1f-%.1f\n', anchors(q), ...
            numel(alignIds{q}), min([alignLen{q}; NaN]), max([alignLen{q}; NaN]));
end
row = zeros(sim.N, 1);  row(sel) = 1:nShow;
r   = row(out.spikes.cell);
plot(out.spikes.t(r > 0), r(r > 0), '.', 'MarkerSize', 5, 'Color', 'k');
xline(ana.tPulse,'LineStyle','-','LineWidth',1,'Color',[1 0 0],'HandleVisibility','Off');
set(gca, 'YDir', 'reverse');  ylim([0, nShow + 1]);
% if ana.rasterBlocks
%     gs = unique(grp, 'stable');                          % groups present, in row order
%     for q = 2:numel(gs)                                  % separators between blocks
%         yb = find(grp == gs(q), 1) - 0.5;
%         plot([0 sim.T], [yb yb], '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, ...
%              'HandleVisibility', 'off');
%     end
% end
ylabel('neurons',"FontSize",12,"Interpreter","latex");
set(gca,'XTickLabel',[]);
%title(sprintf('Raster ($w = %g$, window $\\Delta t = %g$, stride $%g$)', sim.w, ana.win, ana.step));
text(ax(1), -0.035, 1.03, '(a)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'Clipping', 'off');

% --- firing occupancy by pathway ---------------------------------------------------
% Each spike is assigned the pathway that produced it: R = recurrence C_k -> F,
% C = activation Q -> F after a canonical return C_n -> Q, S = activation after a
% sequestration return Y_m -> Q. Components are spikes per neuron (of the whole
% population) per unit time, so they add up to the total (apart from first
% activations from the initial Q).
ax(2) = nexttile; hold on; box on;
pathCols = [1.00,0.30,0.70; 0.3,0.8,0.3; 0.00,0.4,0.70];
pathLeg  = {'Immediate-return', 'Canonical', 'Sequestration'};
for g = 1:3
    plot(S.tc, S.occPath(g, :) / par.gamma(1), '-', 'Color', pathCols(g, :), 'LineWidth', 1.5, ...
         'DisplayName', pathLeg{g});
end
plot(S.tc, S.occupancy / par.gamma(1), ':', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Total');
xline(ana.tPulse,'LineStyle','-','LineWidth',1,'Color',[1 0 0],'HandleVisibility','Off');
ylabel("$p_F$","FontSize",12,"Interpreter","latex");
legend('Location', 'northwest', 'Orientation', 'vertical','FontSize',11);
set(gca,'XTickLabel',[],'YGrid','On');
%title('Firing occupancy $=$ (number of spikes in window)$\,/\,(N\,\Delta t)$, by pathway');
text(ax(2), -0.035, 1.11, '(b)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'Clipping', 'off');

% --- ensemble Fano factor ----------------------------------------------------------
ax(3) = nexttile; hold on; box on;
plot(S.tc, S.fano, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
xline(ana.tPulse,'LineStyle','-','LineWidth',1,'Color',[1 0 0],'HandleVisibility','Off');
ylabel('Fano factor','FontSize',11);
set(gca,'XTickLabel',[],'YGrid','On','YLim',[0.5 1.75],'YTick',[0.5 1 1.5]);

% title('Ensemble Fano factor $\sigma^2 / \mu$ of spike counts in window');
text(ax(3), -0.035, 1.11, '(c)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'Clipping', 'off');

% --- active neurons ----------------------------------------------------------------
ax(4) = nexttile; box on;
plot(S.tc, S.active / sim.N, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
xline(ana.tPulse,'LineStyle','-','LineWidth',1,'Color',[1 0 0],'HandleVisibility','Off');
ylabel("$n_\Delta$","FontSize",12,"Interpreter","latex");
set(gca,'XTickLabel',[],'YGrid','On');

% title('Active proportion $=$ (neurons with $\geq 1$ spike in window)$\,/\,N$');
text(ax(4), -0.035, 1.11, '(d)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'Clipping', 'off');

% --- neurons changing firing rate ----------------------------------------------------
RC = rateChange(S, out.spikes, sim.N, ana.win, ana.rateTol, ana.rateRef, ana.rateSkipSilent);
post = RC.t > ana.tPulse;                   % comparison with the baseline is meaningless before the pulse
fUp = RC.up / sim.N;    fUp(~post) = NaN;
fDn = RC.down / sim.N;  fDn(~post) = NaN;
fSm = RC.same / sim.N;  fSm(~post) = NaN;
ax(5) = nexttile; hold on; box on;
x100 = 100 * ana.rateTol;
plot(RC.t, fUp, '--', 'Color', [1 0.5 0],   'LineWidth', 1.5, ...
     'DisplayName', sprintf('Increasing ($> +%g\\%%$)', x100));
plot(RC.t, fDn, '-', 'Color', [0.6 0.2 1],  'LineWidth', 1.5, ...
     'DisplayName', sprintf('Decreasing ($< -%g\\%%$)', x100));
plot(RC.t, fSm, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 2, ...
     'DisplayName', sprintf('Same (within $\\pm %g\\%%$)', x100));
xline(ana.tPulse,'LineStyle','-','LineWidth',1,'Color',[1 0 0],'HandleVisibility','Off');
ylabel({'fraction'},'FontSize',11);  
set(gca,'XTickLabel',[],'YGrid','On');

l5 = legend('Location', 'northwest', 'Orientation', 'vertical','FontSize',11);
l5.Title.String = {'Fraction of neurons','changing firing rate'};
text(ax(5), -0.035, 1.11, '(e)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'Clipping', 'off');

% Occupancies - Gillespie vs homogeneous mean field
ax(6) = nexttile;
hold on; box on;
cols  = [0 0.5 0; 1 0 0; 0 0 0; 0.00,0.4,0.70];
G = group(out.frac);
for q = 1:4
    plot(t, G(q,:), '-', 'Color', cols(q,:), 'LineWidth', 1, 'DisplayName', names{q});
end
if isfield(out, 'fracHMF')                           % homogeneous mean field (mean parameters)
    GH = group(out.fracHMF);
    for q = 1:4
        plot(t, GH(q,:), ':', 'Color', cols(q,:), 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
xline(ana.tPulse,'LineStyle','-','LineWidth',1,'Color',[1 0 0],'HandleVisibility','Off');
xlabel("Time $t$","FontSize",12,"Interpreter","latex"); 
ylabel("$p_s^*$","FontSize",12,"Interpreter","latex");
legend('Location', 'best', 'Orientation', 'vertical',"FontSize",12,"Interpreter","latex");
set(gca,'XGrid','On','YGrid','On',"FontSize",12,'TickLabelInterpreter','latex');
% set(gca,'YScale','log','YLim',[0.001 1.5],'YTick',[0.001 0.01 0.1 1]);

legend(names, 'Location', 'east', 'Orientation', 'vertical', ...
       'FontSize', 12, 'Interpreter', 'latex');
text(ax(6), -0.035, 1.11, '(f)', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'Clipping', 'off');

linkaxes(ax, 'x');  
set(ax,'FontSize',11,'XGrid','On','XLim',[0 sim.T],'XTick',0:25:sim.T);

% Saving
if save_option
    output_directory = 'figures';

    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end
    

    savefig(fig2, fullfile(output_directory, ...
        "raster_three_pathway.fig"));

    exportgraphics(fig2, fullfile(output_directory, ...
        "raster_three_pathway.png"), ...
        'Resolution', 300);

    exportgraphics(fig2, fullfile(output_directory, ...
        "raster_three_pathway.pdf"), ...
        'ContentType', 'image','Resolution', 300);

    disp('Done');
end


%% =========================================================================
%                              LOCAL FUNCTIONS
%  =========================================================================
function out = simulateCoupledQFCY(par, het, sim)
% Joint SSA of N coupled cells, z = w*N_F/N + u(t), exact via thinning.

N = sim.N;  w = sim.w;  uf = asFunction(sim.u);
[E, exits, labels, iF, nS, P] = buildModel(par, het, N);

% bound of every cell in every state = sum of prefactors of its exits
Bmat = zeros(N, nS);
for s = 1:nS, Bmat(:, s) = sum(E.c(:, exits{s}), 2); end

% transposed parameter tables: cell i's column is contiguous in memory
cT = E.c';  thT = E.th';  muT = E.mu';  isSig = E.isSig;  to = E.to;  typ = E.type;

% initial condition
s = sim.s0(:);  if isscalar(s), s = repmat(s, N, 1); end
nF = sum(s == iF);
b  = Bmat(sub2ind([N nS], (1:N)', s));

% Fenwick tree over the per-cell bounds b
lb   = (1:N)' - bitand((1:N)', (0:N-1)');          % lowest set bit of i
pows = 2.^(floor(log2(N)):-1:0);
tree = fenwickBuild(b, lb);
Btot = sum(b);

% output storage
tgrid = 0:sim.dt:sim.T;  nG = numel(tgrid);
X   = zeros(N, nG, 'uint16');
cnt = zeros(N, 7);
nTr = min(sim.nTrace, N);
trT = cell(nTr, 1);  trS = cell(nTr, 1);
cap  = 8 * N;                                 % spike storage, grown by doubling
spkI = zeros(cap, 1);  spkT = zeros(cap, 1);  spkK = zeros(cap, 1);  nSpk = 0;
spkP = zeros(cap, 1);                         % pathway that produced the spike
lastRet = zeros(N, 1);                        % route of the last return to Q: 2 C, 3 S, 0 none yet
for i = 1:nTr, trT{i} = 0;  trS{i} = s(i); end

t = 0;  g = 1;  nCand = 0;  nAcc = 0;
while Btot > 0
    t = t - log(rand) / Btot;                 % candidate event time
    if t > sim.T, break; end
    nCand = nCand + 1;

    % --- pick a cell proportionally to its bound (Fenwick descent) ------
    target = rand * Btot;  pos = 0;
    for pw = pows
        nxt = pos + pw;
        if nxt <= N && tree(nxt) <= target
            pos = nxt;  target = target - tree(nxt);
        end
    end
    i = min(pos + 1, N);                      % target is now uniform on [0, b(i))

    % --- thinning: accept an exit with probability true rate / bound ------
    z   = w * nF / N + uf(t);                 % current input
    idx = exits{s(i)};
    acc = 0;  cum = 0;
    for q = 1:numel(idx)
        e  = idx(q);
        re = cT(e, i);
        if isSig(e), re = re / (1 + exp(-thT(e, i) * (z - muT(e, i)))); end
        cum = cum + re;
        if target < cum, acc = e; break; end
    end
    if acc == 0, continue; end                % rejected candidate

    % --- accepted: record, jump, update feedback and bounds ---------------
    while g <= nG && tgrid(g) < t
        X(:, g) = s;  g = g + 1;
    end
    sOld = s(i);  sNew = to(acc);  s(i) = sNew;
    nF = nF + (sNew == iF) - (sOld == iF);
    cnt(i, typ(acc)) = cnt(i, typ(acc)) + 1;
    if typ(acc) == 3, lastRet(i) = 2; elseif typ(acc) == 7, lastRet(i) = 3; end
    if i <= nTr, trT{i}(end+1) = t;  trS{i}(end+1) = sNew; end
    if sNew == iF                             % spike = entry into F
        nSpk = nSpk + 1;
        if nSpk > cap
            cap = 2 * cap;  spkI(cap) = 0;  spkT(cap) = 0;  spkK(cap) = 0;  spkP(cap) = 0;
        end
        spkI(nSpk) = i;  spkT(nSpk) = t;  spkK(nSpk) = typ(acc);  % 1 act., 4 recur.
        if typ(acc) == 4, spkP(nSpk) = 1;         % R: recurrence C_k -> F
        else,             spkP(nSpk) = lastRet(i); % activation after C (Cn->Q) or S (Ym->Q) return
        end
    end

    delta = Bmat(i, sNew) - b(i);  b(i) = Bmat(i, sNew);
    k = i;
    while k <= N, tree(k) = tree(k) + delta;  k = k + lb(k); end
    Btot = Btot + delta;

    nAcc = nAcc + 1;
    if mod(nAcc, N) == 0                      % clear floating-point drift
        tree = fenwickBuild(b, lb);  Btot = sum(b);
    end
end
X(:, g:end) = repmat(s, 1, nG - g + 1);

frac = zeros(nS, nG);
for q = 1:nS, frac(q, :) = mean(X == q, 1); end

out.t = tgrid;  out.X = X;  out.frac = frac;  out.counts = cnt;  out.P = P;
out.u = arrayfun(uf, tgrid);
out.z = w * frac(iF, :) + out.u;
out.labels = labels;  out.edges = E;
out.traces = struct('t', cellfun(@(v) [v, sim.T], trT, 'UniformOutput', false), ...
                    's', cellfun(@(v) [v, v(end)], trS, 'UniformOutput', false));
out.nCandidates = nCand;  out.nAccepted = nAcc;
out.spikes.cell = spkI(1:nSpk);  out.spikes.t = spkT(1:nSpk);  out.spikes.type = spkK(1:nSpk);
out.spikes.path = spkP(1:nSpk);   % 1 R, 2 C, 3 S, 0 activation before any return (initial Q)
out.typeNames = {'activation Q->F', 'canonical step', 'canonical return Cn->Q', ...
                 'recurrence Ck->F', 'sequestration entry Cl->Y1', ...
                 'sequestration step', 'sequestration return Ym->Q'};

mf = 'none';
if isfield(sim, 'meanField'), mf = sim.meanField; end
if islogical(mf) || isnumeric(mf)                    % backward compatible: true -> 'closed'
    if mf, mf = 'closed'; else, mf = 'none'; end
end
mf = lower(char(mf));
s0 = sim.s0(:);  if isscalar(s0), s0 = repmat(s0, N, 1); end
if any(strcmp(mf, {'closed', 'both'}))               % self-consistent N -> infinity limit
    [out.fracMF, out.zMF] = meanFieldQuenched(E, s0, nS, iF, w, uf, tgrid, []);
end
if any(strcmp(mf, {'driven', 'both'}))               % driven by the SSA's realised p_F(t)
    [out.fracDMF, out.zDMF] = meanFieldQuenched(E, s0, nS, iF, w, uf, tgrid, frac(iF, :));
end
if isfield(sim, 'meanFieldHom') && sim.meanFieldHom  % homogeneous: mean parameters, one cell
    rngState = rng;                                  % keep the random stream unchanged
    het0 = structfun(@(x) 0 * x, het, 'UniformOutput', false);
    E0   = buildModel(par, het0, 1);                 % all CVs/SDs = 0 -> exactly the mean values
    rng(rngState);
    p0   = accumarray(s0, 1, [nS 1])' / N;           % initial state distribution of the population
    [out.fracHMF, out.zHMF] = meanFieldQuenched(E0, 1, nS, iF, w, uf, tgrid, [], p0);
end
end

function [pbar, zMF] = meanFieldQuenched(E, s0, nS, iF, w, uf, tgrid, pFdrive, p0)
% Mean-field master equations for the sampled parameters:
%   dp_i/dt = p_i * G_i(z)
%   pFdrive empty: z = w * mean_i p_i,F + u(t)        (closed loop, N -> infinity)
%   pFdrive given: z = w * pFdrive(t) + u(t)          (driven, pFdrive on tgrid,
%                                                      linearly interpolated)
% p0 (optional, N x nS): initial probabilities instead of one-hot rows at s0.
% Integrated with adaptive ode45 over chunks of the output grid (low memory).
N  = size(E.c, 1);  nE = numel(E.from);  nG = numel(tgrid);
D  = zeros(nE, nS);                          % edge -> state incidence
for e = 1:nE, D(e, E.from(e)) = -1;  D(e, E.to(e)) = D(e, E.to(e)) + 1; end
sg = E.isSig;
if isempty(pFdrive)
    zfun = @(t, p) w * mean(p(:, iF)) + uf(t);
else
    dtg  = tgrid(2) - tgrid(1);
    zfun = @(t, p) w * gridInterp(pFdrive, t, dtg) + uf(t);
end
f  = @(t, y) reshape(mfRHS(t, reshape(y, N, nS), E.c, E.th, E.mu, sg, ...
                           E.from, D, zfun), [], 1);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', tgrid(2) - tgrid(1));

if nargin >= 9 && ~isempty(p0)
    p = p0;
else
    p = zeros(N, nS);  p(sub2ind([N nS], (1:N)', s0)) = 1;
end
pbar = zeros(nS, nG);  pbar(:, 1) = mean(p, 1)';
y = p(:);  g0 = 1;  chunk = 40;
while g0 < nG
    g1 = min(g0 + chunk, nG);
    [~, Y] = ode45(f, tgrid(g0:g1), y, opts);
    if g1 - g0 == 1, Y = Y([1 end], :); end  % 2-point tspan returns all steps
    for r = 2:size(Y, 1)
        pbar(:, g0 + r - 1) = mean(reshape(Y(r, :), N, nS), 1)';
    end
    y = Y(end, :)';  g0 = g1;
end
zMF = w * pbar(iF, :) + arrayfun(uf, tgrid);
end

function dp = mfRHS(t, p, c, th, mu, sg, from, D, zfun)
z = zfun(t, p);
r = c;
r(:, sg) = c(:, sg) ./ (1 + exp(-th(:, sg) .* (z - mu(:, sg))));
dp = (p(:, from) .* r) * D;
end

function v = gridInterp(x, t, dt)
% linear interpolation of x (sampled at 0, dt, 2*dt, ...) at time t
k = min(max(floor(t / dt) + 1, 1), numel(x) - 1);
a = min(max(t / dt - (k - 1), 0), 1);
v = (1 - a) * x(k) + a * x(k + 1);
end

function R = rateChange(S, spk, N, win, tol, ref, skipSilent)
% Neurons firing more, fewer or about the same number of spikes in each window
% as expected from their own baseline rate.
%   r0       = spikes in the baseline [ref(1), ref(2)) / its length  (can be zero)
%   expected = r0 * win   (expected spikes in a window of length win)
%   increasing : count > (1 + tol) * expected
%   decreasing : count < (1 - tol) * expected
%   same       : otherwise
% skipSilent: neurons with zero count and zero baseline are left out of 'same'.
r0       = accumarray(spk.cell(spk.t >= ref(1) & spk.t < ref(2)), 1, [N 1]) / diff(ref);
expected = r0 * win;                    % N x 1
W        = S.counts;                    % N x nWin spike counts per window
up   = W > (1 + tol) * expected;
down = W < (1 - tol) * expected;
same = ~up & ~down;
if skipSilent, same = same & (W > 0 | expected > 0); end
R.t = S.tc;
R.up = sum(up, 1);  R.down = sum(down, 1);  R.same = sum(same, 1);
R.refText = sprintf('baseline $[%g, %g)$', ref(1), ref(2));
end

function S = spikeStatistics(spk, N, T, win, step, minISI)
% Windowed ensemble statistics of the spike trains (spike = entry into F).
% Windows [t0, t0 + win) with t0 = 0, step, 2*step, ...; S.tc = window ENDS, so each
% point only uses spikes up to that time (features line up with their cause).
%   S.counts     N x nWin spike counts per neuron and window
%   S.total      total spikes in each window
%   S.occupancy  firing occupancy = total / (N * win)  (spikes per neuron per unit time)
%   S.occPath    3 x nWin occupancy split by spike pathway (R, C, S), same units
%   S.active     number of neurons with >= 1 spike in the window
%   S.fano       Fano factor var/mean of the per-neuron counts across the ensemble
%   S.cvISI      CV of all inter-spike intervals (pooled over neurons) that
%                end inside the window; NaN if fewer than minISI intervals
k = round(win / step);
assert(k >= 1 && abs(k*step - win) < 1e-9, 'ana.win must be an integer multiple of ana.step.');
nb = floor(T / step + 1e-9);                         % fine bins of width step
assert(nb >= k, 'Window longer than the simulation.');

b  = floor(spk.t / step) + 1;  kp = b <= nb;
Cb = accumarray([spk.cell(kp), b(kp)], 1, [N nb]);   % counts per neuron per fine bin
W  = winSum(Cb, k);                                  % counts per neuron per window
nw = size(W, 2);

S.tc        = (0:nw-1) * step + win;                % trailing window: plotted at its end
S.counts    = W;
S.total     = sum(W, 1);
S.occupancy = S.total / (N * win);
S.active    = sum(W > 0, 1);
mu = mean(W, 1);  vr = var(W, 0, 1);
S.fano      = vr ./ mu;   S.fano(mu == 0) = NaN;

% pooled inter-spike intervals, assigned to the window containing their end
[~, ord] = sortrows([spk.cell, spk.t]);
c = spk.cell(ord);  ts = spk.t(ord);
same = [false; c(2:end) == c(1:end-1)];
isi  = [NaN; diff(ts)];   isi = isi(same);   bi = floor(ts(same) / step) + 1;
kp   = bi <= nb;  isi = isi(kp);  bi = bi(kp);
n  = winSum(accumarray(bi, 1,      [nb 1])', k);
s1 = winSum(accumarray(bi, isi,    [nb 1])', k);
s2 = winSum(accumarray(bi, isi.^2, [nb 1])', k);
m1 = s1 ./ n;
v  = (s2 - n .* m1.^2) ./ (n - 1);
S.cvISI = sqrt(max(v, 0)) ./ m1;   S.cvISI(n < max(minISI, 2)) = NaN;
S.nISI  = n;

% occupancy split by the pathway that produced each spike (spk.path: 1 R, 2 C, 3 S);
% rows sum to S.occupancy except for spikes with path 0 (activation from the initial Q)
S.occPath = zeros(3, nw);
if isfield(spk, 'path')
    pth = spk.path(:);  kb = b <= nb;
    for g = 1:3
        sel = kb & pth == g;
        S.occPath(g, :) = winSum(accumarray(b(sel), 1, [nb 1])', k) / (N * win);
    end
end
end

function cls = pathwayClass(cnt, ratio)
% Pathway class of each cell from its event counts over [0, T]
% (columns of cnt: 3 canonical return Cn->Q, 4 recurrence Ck->F, 5 sequestration Cl->Y1)
%   3 S  sequestered at least once (rare, so it takes priority)
%   1 R  otherwise, recurrences >= ratio * canonical returns (and at least one recurrence)
%   2 C  otherwise, i.e. recurrences < ratio * canonical returns
%   0    no pathway exit at all during the run
rec = cnt(:, 4);  can = cnt(:, 3);  seq = cnt(:, 5) > 0;
cls = zeros(size(cnt, 1), 1);
isR = ~seq & rec > 0 & rec >= ratio * can;
cls(isR)                         = 1;
cls(~seq & ~isR & can > 0)       = 2;
cls(seq)                         = 3;
end

function [ids, len] = alignedSilences(spk, N, T, anchors, tol, nPer)
% For each anchor time a, the nPer neurons with the longest silent gap that
% starts in [a, a + tol), i.e. last spike in that interval and no further spike
% for as long as possible (gaps still open at T end at T). A neuron is used in
% at most one group. ids{q}, len{q}: neurons (longest gap first), gap lengths.
[~, ord] = sortrows([spk.cell(:), spk.t(:)]);
c  = spk.cell(ord);  ts = spk.t(ord);
nx = [ts(2:end); T];
nx([c(2:end) ~= c(1:end-1); true]) = T;         % last spike of each neuron: gap runs to T
gap  = nx - ts;
used = false(N, 1);
ids = cell(numel(anchors), 1);  len = ids;
for q = 1:numel(anchors)
    m    = ts >= anchors(q) & ts < anchors(q) + tol;
    best = accumarray(c(m), gap(m), [N 1], @max, 0);   % longest such gap per neuron
    best(used) = 0;
    [v, o] = sort(best, 'descend');
    kq = min(nPer, sum(v > 0));
    ids{q} = o(1:kq);  len{q} = v(1:kq);  used(ids{q}) = true;
end
end

function [sel, grp] = selectRasterCells(cls, tot, nShow, frac, blocks)
% Random subset of nShow cells stratified by pathway class with target shares
% frac = [R C S]. If a class has too few cells, all of them are taken and the
% shortfall is filled at random from the remaining classified cells (then from
% unclassified ones). Rows are sorted by spike count, within each class block
% if blocks is true. grp = row class (1 R, 2 C, 3 S, 4 other).
frac  = frac(:)' / sum(frac);
quota = frac * nShow;  nq = floor(quota);                 % largest-remainder rounding
[~, ix] = sort(quota - nq, 'descend');  d = nShow - sum(nq);
nq(ix(1:d)) = nq(ix(1:d)) + 1;

sel = zeros(0, 1);
for g = 1:3
    pool = find(cls == g);  kq = min(nq(g), numel(pool));
    sel  = [sel; pool(randperm(numel(pool), kq))];        %#ok<AGROW>
end
short = nShow - numel(sel);
if short > 0
    rest = setdiff(find(cls > 0), sel);
    if numel(rest) < short, rest = setdiff((1:numel(cls))', sel); end
    sel  = [sel; rest(randperm(numel(rest), min(short, numel(rest))))];
    warning('Raster: not enough cells in some pathway class(es); %d rows filled from other classes.', short);
end

grp = cls(sel);  grp(grp == 0) = 4;
if blocks, [~, o] = sortrows([grp, -tot(sel)]);
else,      [~, o] = sort(tot(sel), 'descend');
end
sel = sel(o);  grp = grp(o);
end

function y = winSum(x, k)
% sums over k consecutive columns (sliding window) via cumulative sums
c = cumsum([zeros(size(x, 1), 1), x], 2);
y = c(:, k+1:end) - c(:, 1:end-k);
end

function [E, exits, labels, iF, nS, P] = buildModel(par, het, N)
n = par.n;  m = par.m;  kk = par.k(:)';  ll = par.l(:)';  nk = numel(kk);  nl = numel(ll);
assert(m >= 1, 'Need at least one sequestration stage (m >= 1).');
assert(numel(par.gamma) == n + 1, 'par.gamma must have n+1 entries.');
assert(numel(par.alpha) == m,     'par.alpha must have m entries.');
assert(all(kk >= 1 & kk <= n) && all(ll >= 1 & ll <= n), 'k, l must be in 1..n.');
par.beta  = expandTo(par.beta,  nk, 'beta');   par.kappa = expandTo(par.kappa, nl, 'kappa');
par.thP   = expandTo(par.thP,   nk, 'thP');    par.thD   = expandTo(par.thD,   nl, 'thD');
par.muP   = expandTo(par.muP,   nk, 'muP');    par.muD   = expandTo(par.muD,   nl, 'muD');

iQ = 1;  iF = 2;  iC = 2 + (1:n);  iY = 2 + n + (1:m);  nS = 2 + n + m;
labels = [{'Q', 'F'}, ...
          arrayfun(@(j) sprintf('C%d', j), 1:n, 'UniformOutput', false), ...
          arrayfun(@(j) sprintf('Y%d', j), 1:m, 'UniformOutput', false)];

P = drawQuenched(par, het, N);               % quenched: drawn once, never redrawn

% type: 1 activation, 2 canonical step, 3 canonical return Cn->Q,
%       4 recurrence Ck->F, 5 sequestration entry Cl->Y1,
%       6 sequestration step, 7 sequestration return Ym->Q
E.from = zeros(1, 0);  E.to = zeros(1, 0);  E.type = zeros(1, 0);  E.isSig = false(1, 0);
E.c = zeros(N, 0);  E.th = zeros(N, 0);  E.mu = zeros(N, 0);
O = ones(N, 1);  Z = zeros(N, 1);

E = addEdge(E, iQ, iF, 1, true, O, P.thA, P.muA);             % A(z), A_max = 1
E = addEdge(E, iF, iC(1), 2, false, P.gamma(:, 1), Z, Z);     % gamma_0
for j = 1:n
    if j < n, E = addEdge(E, iC(j), iC(j+1), 2, false, P.gamma(:, j+1), Z, Z);
    else,     E = addEdge(E, iC(n), iQ,      3, false, P.gamma(:, n+1), Z, Z);
    end
end
for q = 1:nk
    E = addEdge(E, iC(kk(q)), iF, 4, true, P.beta(:, q), P.thP(:, q), P.muP(:, q));
end
for q = 1:nl
    E = addEdge(E, iC(ll(q)), iY(1), 5, true, P.kappa(:, q), P.thD(:, q), P.muD(:, q));
end
for j = 1:m
    if j < m, E = addEdge(E, iY(j), iY(j+1), 6, false, P.alpha(:, j), Z, Z);
    else,     E = addEdge(E, iY(m), iQ,      7, false, P.alpha(:, m), Z, Z);
    end
end
exits = arrayfun(@(s) find(E.from == s), 1:nS, 'UniformOutput', false);
end

function P = drawQuenched(par, het, N)
% log-normal keeps the population mean (and sign) with the requested CV
LN = @(mu, cv) mu .* exp(sqrt(log(1 + cv^2)) * randn(N, numel(mu)) - log(1 + cv^2)/2);
GS = @(mu, sd) mu + sd * randn(N, numel(mu));
P.gamma = LN(par.gamma, het.gamma);   P.alpha = LN(par.alpha, het.alpha);
P.beta  = LN(par.beta,  het.beta);    P.kappa = LN(par.kappa, het.kappa);
P.thA = LN(par.thA, het.thA);   P.muA = GS(par.muA, het.muA);
P.thP = LN(par.thP, het.thP);   P.muP = GS(par.muP, het.muP);
P.thD = LN(par.thD, het.thD);   P.muD = GS(par.muD, het.muD);
end

function E = addEdge(E, from, to, type, isSig, c, th, mu)
E.from(end+1) = from;  E.to(end+1) = to;  E.type(end+1) = type;  E.isSig(end+1) = isSig;
E.c(:, end+1) = c;     E.th(:, end+1) = th;  E.mu(:, end+1) = mu;
end

function tree = fenwickBuild(b, lb)
N = numel(b);  tree = b;
for i = 1:N
    j = i + lb(i);
    if j <= N, tree(j) = tree(j) + tree(i); end
end
end

function f = asFunction(u)
if isa(u, 'function_handle'), f = u; else, u0 = u; f = @(t) u0 + 0*t; end
end

function v = expandTo(v, len, name)
if isscalar(v), v = repmat(v, 1, len); end
assert(numel(v) == len, '%s must be a scalar or have %d entries.', name, len);
v = v(:)';
end

