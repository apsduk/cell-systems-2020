%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all');

%% %%%%% SET UP MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runintmax = 120;
tind = 6;

odesolver = @ode23s;
odeoptions = odeset(); % odeset('AbsTol', 1e-6, 'RelTol', 1e-6);

linewidth = 1.5;
fontsize = 10;
outlinewidth = 1;

% cmap = [0 0 0; lines(7)];
cmap = lines(7);
cmap = [0.5 0.5 0.5; cmap(5,:)];

fplot = figure('Units', 'centimeters', 'Position', [5 5 25 30]);

%% %%%%% SET UP MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model{1}.modelid = 'uber_or';
model{1}.Y0 = zeros(25,1);
model{1}.f = @(T, Y, X, Z) uber_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z(1), Z(2), X(44), X(32)], X(33:43), X(45:55));
model{1}.m1idx =  9;
model{1}.p1idx = 11;

model{2}.modelid = 'frag_or';
model{2}.Y0 = zeros(25,1);
model{2}.f = @(T, Y, X, Z) frag_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [ Z(1), Z(2), X(10), X(19), X(42)], X(32:41), X(43:53));
model{2}.m1idx = 14;
model{2}.p1idx = 16;

model{3}.modelid = 'sqtr_or';
model{3}.Y0 = zeros(23,1);
model{3}.f = @(T, Y, X, Z) sqtr_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z( 1), Z( 2), X(40), X(32)], X(33:39), X(41:51));
model{3}.m1idx =  9;
model{3}.p1idx = 11;

model{4}.modelid = 'op_or';
model{4}.Y0 = zeros(20,1);
model{4}.f = @(T, Y, X, Z) op_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z(1), Z(2), X(32)], [], X(33:43));
model{4}.m1idx =  9;
model{4}.p1idx = 11;

model{5}.modelid = 'var_frag_or';
model{5}.Y0 = zeros(25,1);
model{5}.f = @(T, Y, X, Z) var_frag_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [ Z(1), Z(2), X(10), X(19), X(42)], X(32:41), X(43:53));
model{5}.m1idx = 14;
model{5}.p1idx = 16;

model{6}.modelid = 'ortho_op_or';
model{6}.Y0 = zeros(20,1);
model{6}.f = @(T, Y, X, Z) ortho_op_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z(1), Z(2), X(32)], [], X(33:43));
model{6}.m1idx =  9;
model{6}.p1idx = 11;

%% %%%%% CHOOSE MODELS TO SIMULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
topo(1) = 2; mlim(1,:) = [0.95, 1.05]; plim(1,:) = [0.7, 1.1]; controllername{1}{1} = 'FRAG-OR';
topo(2) = 5; mlim(2,:) = [0.95, 1.05]; plim(2,:) = [0.7, 1.1]; controllername{1}{2} = 'var-FRAG-OR';

%% %%%%% MAKE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fplot.Number);

sp{1} = subplot(5, 2, [3,  5]); box('on'); hold('on'); grid('on'); pbaspect([1 1 1]);
sp{2} = subplot(5, 2, [4,  6]); box('on'); hold('on'); grid('on'); pbaspect([1 1 1]);
sp{3} = subplot(5, 2, [7,  9]); box('on'); hold('on'); pbaspect([1 1 1]);
sp{4} = subplot(5, 2, [8, 10]); box('on'); hold('on'); pbaspect([1 1 1]);
ax{1} = axes('Position', [0.335, 0.1285, 0.12, 0.12]); box('on'); hold('on'); % pbaspect([1 1 1]);
ax{2} = axes('Position', [0.775, 0.1285, 0.12, 0.12]); box('on'); hold('on'); % pbaspect([1 1 1]);

figure(fplot.Number);
for m = 1:length(topo)
    
    % make submodel
    submodel = model{topo(m)}
    
    %%%%%%%% PLOT ALL OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load performance metrics
    load(['model_',submodel.modelid,'_different_metrics.mat'], 'couplingsort', 'coupling2dscore', 'proteinfinaloutput');
    
    % only show top controllers
    couplingsort = couplingsort(1:200);
    coupling2dscore = coupling2dscore(1:200);
    proteinfinaloutput = proteinfinaloutput(1:200);
    scatter(sp{1}, couplingsort, coupling2dscore, 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
    set(sp{1}, 'YScale', 'log', 'XLim', [1, length(couplingsort)], 'XTick', [1, 200]);
    xlabel(sp{1}, 'Design index (sorted)')
    ylabel(sp{1}, '2D score');
    set(sp{1}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    scatter(sp{2}, proteinfinaloutput, coupling2dscore, 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
    xlim([1e3, 1e6]);
    set(sp{2}, 'XScale', 'log', 'YScale', 'log', 'XLim', [1e3, 1e6]);
    xlabel(sp{2}, 'Final protein output (nM)');
    ylabel(sp{2}, '2D score');
    set(sp{2}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    %%%%%%%% SIMULATE PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load parameters
    load(['model_',submodel.modelid,'_best_controllers.mat'], 'bestX0');
    
    % select best parameters
    X = bestX0(1,:);
    
    % make ode functions
    f_0 = @(T, Y) submodel.f(T, Y, X, [   0,    0]);
    f_1 = @(T, Y) submodel.f(T, Y, X, [X(1),    0]);
    f_2 = @(T, Y) submodel.f(T, Y, X, [X(1), X(1)]);
    
    % simulate odes
    [T_0, Y_0] = odesolver(f_0, [- runintmax, 0], submodel.Y0, odeoptions);
    [T_0, Y_0] = odesolver(f_1, [- runintmax, 0], Y_0(end,:)', odeoptions);
    [T_1, Y_1] = odesolver(f_1, [0, tind], Y_0(end,:)', odeoptions);
    [T_2, Y_2] = odesolver(f_2, [max(T_1) + (1e-6), max(T_1) + 3*tind], Y_1(end,:)', odeoptions);
    T = [T_0; T_1; T_2];
    Y = [Y_0; Y_1; Y_2];
    
    plot(sp{3}, T, Y(:,submodel.m1idx)./Y_1(end,submodel.m1idx), '-', 'Color', cmap(m,:), 'LineWidth', 2);
    set(sp{3}, 'YLim', [0.950, 1.01], 'XLim', [0, max(T)], 'XTick', [0:6:max(T)]);
    ylabel(sp{3}, '\Delta TX');
    xlabel(sp{3}, 'Time (h)');
    set(sp{3}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    bar(ax{1}, m, Y(end, submodel.m1idx), 'FaceColor', cmap(m,:));
    set(ax{1}, 'YScale', 'log', 'YLim', [1e3 10^6.5], 'YTick', [1e3, 1e4, 1e5, 1e6], 'XTick', [], 'XLim', [0.5, 2.5]);
    ylabel(ax{1}, 'mRNA output');
    set(ax{1}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    plot(sp{4}, T, Y(:,submodel.p1idx)./Y_1(end,submodel.p1idx), '-', 'Color', cmap(m,:), 'LineWidth', 2);
    set(sp{4}, 'YLim', [0.7 1.1], 'XLim', [0, max(T)], 'XTick', [0:6:max(T)]);
    ylabel(sp{4}, '\Delta TL');
    xlabel(sp{4}, 'Time (h)');
    set(sp{4}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    bar(ax{2}, m, Y(end, submodel.p1idx), 'FaceColor', cmap(m,:));
    set(ax{2}, 'YScale', 'log', 'YLim', [1e3 1e6], 'YTick', [1e3, 1e4, 1e5, 1e6], 'XTick', [], 'XLim', [0.5, 2.5]);
    ylabel(ax{2}, 'Protein output');
    set(ax{2}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    legend(sp{1}, controllername{1}, 'NumColumns', 3, 'FontSize', 14);
    
end
