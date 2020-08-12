%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all'); addpath('libmodels');

hname = pwd;

%% %%%%% SET UP MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runintmax = 120;
tind = 12;
tmax = 12;

odesolver = @ode23s;
odeoptions = odeset(); % odeset('AbsTol', 1e-6, 'RelTol', 1e-6);

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

%% %%%%% PLOT TOPOLOGY SCORES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load performance metrics
load(['model_ortho_op_or_different_metrics.mat'], 'couplingsort', 'coupling2dscore', 'proteinfinaloutput');
couplingsort = couplingsort(1:200);
coupling2dscore = coupling2dscore(1:200);
proteinfinaloutput = proteinfinaloutput(1:200);

fplot = figure('Units', 'centimeters', 'Position', [5 5 20 10]);
subplot(1, 2, 1); hold('on'); grid('on'); box('on'); pbaspect([1 1 1]);
scatter(couplingsort, coupling2dscore, 'o', 'MarkerFaceColor', lines(1), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
set(gca, 'YScale', 'log', 'XLim', [1, length(couplingsort)], 'XTick', [1, length(couplingsort)], 'LineWidth', 1, 'FontWeight', 'bold');
xlabel('Design index (sorted)');
ylabel('2D score');
subplot(1, 2, 2); hold('on'); grid('on'); box('on'); pbaspect([1 1 1]);
scatter(proteinfinaloutput, coupling2dscore, 'o', 'MarkerFaceColor', lines(1), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [1e3, 1e6], 'LineWidth', 1, 'FontWeight', 'bold');
xlabel('Final protein output (nM)');
ylabel('2D score');

%% %%%%% SIMULATE DYNAMICS OF ORTHO-OP-OR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select model
submodel = model{6};

% load parameters
load(['model_ortho_op_or_best_controllers.mat'], 'bestX0');

% select best parameters
X = bestX0(1,:);

% make ode functions
f_0 = @(T, Y) submodel.f(T, Y, X, [   0,    0]);
f_1 = @(T, Y) submodel.f(T, Y, X, [X(1),    0]);
f_2 = @(T, Y) submodel.f(T, Y, X, [X(1), X(1)]);

% simulate odes
[T_0, Y_0] = odesolver(f_0, [- runintmax, 0], submodel.Y0, odeoptions);
[T_0, Y_0] = odesolver(f_0, [- tind, 0.5*tind], Y_0(end,:)', odeoptions);
[T_1, Y_1] = odesolver(f_1, [0.5*tind + (1e-6), max(T_0) + tind], Y_0(end,:)', odeoptions);
[T_2, Y_2] = odesolver(f_2, [max(T_1) + (1e-6), max(T_1) + tmax], Y_1(end,:)', odeoptions);
T = [T_0; T_1; T_2];
Y = [Y_0; Y_1; Y_2];

% model dynamics
xP = Y(:, 1); mP = Y(:, 2); cP = Y(:, 3); oP = Y(:, 4);
xR = Y(:, 5); rR = Y(:, 6); oR = Y(:, 7); x1 = Y(:, 8);
m1 = Y(:, 9); c1 = Y(:,10); p1 = Y(:,11); x2 = Y(:,12);
m2 = Y(:,13); c2 = Y(:,14); p2 = Y(:,15);
kR = Y(:,16); xF = Y(:,17); mF = Y(:,18); cF = Y(:,19); pF = Y(:,20);

%% %%%%% PLOT OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = lines(4);

figure('Units', 'centimeters', 'Position', [5 5 25 25]);
subplot(2, 2, 1); box('on'); hold('on');
plot(T - 6, m1, '-', 'Color', cmap(1,:), 'LineWidth', 1.5);
plot(T - 6, m2, '-', 'Color', cmap(2,:), 'LineWidth', 1.5);
ylim([0, 1.1*max(m1)]); pbaspect([1 1 1]);
xlim([-6, max(T)-6]);
xticks([0:6:max(T)-6]);
ylabel('mRNA output (nM)');
xlabel('Time (h)');
legend('m_1','m_2','Location', 'southeast');
set(gca, 'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');
subplot(2, 2, 2); box('on'); hold('on');
plot(T - 6, p1, '-', 'Color', cmap(1,:), 'LineWidth', 1.5);
plot(T - 6, p2, '-', 'Color', cmap(2,:), 'LineWidth', 1.5);
ylim([0, 1.1*max(p1)]); pbaspect([1 1 1]);
xlim([-6, max(T)-6]);
xticks([0:6:max(T)-6]);
ylabel('Protein output (nM)');
xlabel('Time (h)');
legend('p_1','p_2','Location', 'southeast');
set(gca, 'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');

idx = sum(T < T_0(end));

subplot(2, 2, 3); box('on'); hold('on');
plot(T - 6, oP./oP(idx), '-', 'Color', cmap(4,:), 'LineWidth', 1.5);
plot(T - 6, (oP + xP + xR + x1 + x2)./(oP(idx) + xP(idx) + xR(idx) + x1(idx) + x2(idx)), '-', 'Color', 'k', 'LineWidth', 1.5);
xlim([-6, max(T)-6]);
xticks([0:6:max(T)-6]);
ylim([0.9, 1.1]);
yticks([0.9:0.05:1.1]);
ylabel('Transcriptional activity');
xlabel('Time (h)');
legend('Free o-RNAP','Total o-RNAP','Location', 'northwest');
set(gca, 'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');

subplot(2, 2, 4); box('on'); hold('on');
plot(T - 6, oR./oR(idx), '-', 'Color', cmap(3,:), 'LineWidth', 1.5);
plot(T - 6, (oR + cP + cF + c1 + c2)./(oR(idx) + cP(idx) + cF(idx) + c1(idx) + c2(idx)), '-', 'Color', 'k', 'LineWidth', 1.5);
xlim([-6, max(T)-6]);
xticks([0:6:max(T)-6]);
ylim([0.5, 3.6]);
yticks([0.5:0.5:3.5]);
ylabel('Translational activity');
xlabel('Time (h)');
legend('Free o-ribosomes','Total o-ribosomes','Location', 'northwest');
set(gca, 'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');
