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

cmap = lines(2);

fplot = figure('Units', 'centimeters', 'Position', [0 0 24 32]);

%% %%%%% LOAD PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% base parameters
oleak = 1;
xif1 = 30; xir1 = 6000; tau1 = 250; deltam1 = 20; betaf1 = 1e2; betar1 = 1e6; gamma1 = 300; deltap1 = 0;
xif2 = 30; xir2 = 6000; tau2 = 250; deltam2 = 20; betaf2 = 1e2; betar2 = 1e6; gamma2 = 300; deltap2 = 0;
xifP = 30; xirP = 6000; tauP = 250; deltamP = 20; betafP = 1e2; betarP = 1e6; gammaP = 300; deltaoP = 0;
xifR = 30; xirR = 6000; tauR = 190; deltarR = 20; rhof = 0.9; rhor = 24.8;
xifF = 30; xirF = 6000; tauF = 250; deltamF = 20; betafF = 1e2; betarF = 1e6; gammaF = 300; deltapF = 0; alphafF = 1; alpharF = 0.02; etaF = 4;
xifQ = 30; xirQ = 6000; tauQ = 250; deltamQ = 20; betafQ = 1e2; betarQ = 1e6; gammaQ = 300; deltapQ = 0; alphafQ = 1; alpharQ = 167; etaQ = 1;
lambda = 1; alphaP = 250; alphaR = 2500; phiP = 10; phiR = 20;
omegagP = 100; omegagR = 100; omegag1 = 500; omegag2 = 500; omegagF = 10; omegagQ = 10;

% ornap dissociation constant
kX0 = 200;

% calculate xif1 and xif2
xif1 = (xir1 + phiR*tau1)./kX0;
xif2 = xif1;

% make circuit vectors
% betaf2 = 0;
circPR = [xif1 xir1 tau1 deltam1 betaf1 betar1 gamma1 deltap1 xif2 xir2 tau2 deltam2 betaf2 betar2 gamma2 deltap2];

% host prameters with o leak
hostPR = [lambda alphaP alphaR phiP phiR oleak];

%% %%%%% ISOLATED CONTROLLERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% TRANSCRIPTIONAL CONTROLLER 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specific circuit parameters
omegagP = 100;
omegagQ = 10;
betafP = 100.03;
betafQ = 100.03;
xifP = 55;
xifQ = 55;
alpharQ = 167;
etaQ = 1;

% make vectors
ornapPR1 = [xifP xirP tauP deltamP betafP betarP gammaP deltaoP];
rnapcontPR1 = [xifQ xirQ tauQ deltamQ betafQ betarQ gammaQ deltapQ alphafQ alpharQ etaQ];

%%%%%%%% TRANSCRIPTIONAL CONTROLLER 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xifA = xifP; xirA = xirP; tauA = tauP; deltamA = deltamP; betafA = betafP; betarA = betarP; gammaA = gammaP; deltapA = 0;
xifC = xifP; xirC = xirP; tauC = tauP; deltamC = deltamP; betafC = betafP; betarC = betarP; gammaC = gammaP;

% specific circuit parameters
xifC = 30;
xifA = 30;
betafC = 100.03;
betafA = 100.03;
alphafA = 1.62;  %  :: From Segall-Shapiro2014 ka = 4.5e5 M/s
alpharA = 0.72;  %  :: From Segall-Shapiro2014 kd = 2.0e-4 1/s
omegagC = 500;

% make vectors
ornapPR2 = [xifC xirC tauC deltamC betafC betarC gammaC deltaoP];
rnapcontPR2 = [xifA xirA tauA deltamA betafA betarA gammaA deltapA alphafA alpharA];

%%%%%%%% TRANSLATIONAL CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specific controller parameters
omegagR = 100;
omegagF = 10;
xifR = 20;
alpharF = 0.02;
etaF = 4;
xifF = 30;
betafF = 100.03;

% make vectors
oriboPR = [xifR xirR tauR deltarR rhof rhor];
ribocontPR = [xifF xirF tauF deltamF betafF betarF gammaF deltapF alphafF alpharF etaF];

%%%%%%%% MAKE ODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_1_0 = @(T,Y) uber_or(T, Y, circPR, ornapPR1, oriboPR, hostPR, [omegagP omegagR       0       0 omegagF omegagQ], rnapcontPR1, ribocontPR);
f_1_1 = @(T,Y) uber_or(T, Y, circPR, ornapPR1, oriboPR, hostPR, [omegagP omegagR omegag1       0 omegagF omegagQ], rnapcontPR1, ribocontPR);
f_1_2 = @(T,Y) uber_or(T, Y, circPR, ornapPR1, oriboPR, hostPR, [omegagP omegagR omegag1 omegag2 omegagF omegagQ], rnapcontPR1, ribocontPR);

f_2_0 = @(T,Y) frag_or(T, Y, circPR, ornapPR2, oriboPR, hostPR, [      0       0 omegagC omegagR omegagF], rnapcontPR2, ribocontPR);
f_2_1 = @(T,Y) frag_or(T, Y, circPR, ornapPR2, oriboPR, hostPR, [omegag1       0 omegagC omegagR omegagF], rnapcontPR2, ribocontPR);
f_2_2 = @(T,Y) frag_or(T, Y, circPR, ornapPR2, oriboPR, hostPR, [omegag1 omegag2 omegagC omegagR omegagF], rnapcontPR2, ribocontPR);

%%%%%%%% SIMULATE FIRST COMBINED CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
Y0 = zeros(25,1);
[T_0, Y_0] = odesolver(f_1_0, [- runintmax, 0], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_1_1, [- runintmax, 0], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_1_1, [0, 2*tind], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_1_2, [max(T_1) + (1e-6), max(T_1) + 6*tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

xP = Y(:, 1); mP = Y(:, 2); cP = Y(:, 3); oP = Y(:, 4);
xR = Y(:, 5); rR = Y(:, 6); oR = Y(:, 7);
x1 = Y(:, 8); m1 = Y(:, 9); c1 = Y(:,10); p1 = Y(:,11);
x2 = Y(:,12); m2 = Y(:,13); c2 = Y(:,14); p2 = Y(:,15);
kR = Y(:,16); xF = Y(:,17); mF = Y(:,18); cF = Y(:,19); pF = Y(:,20);
kP = Y(:,21); xQ = Y(:,22); mQ = Y(:,23); cQ = Y(:,24); pQ = Y(:,25);

figure(fplot.Number);
ax{1} = axes('Position', [0.70, 0.75, 0.2, 0.15]); box('on'); hold('on'); % pbaspect([1 1 1]);
yyaxis('left');
plot(T, m1, ':', 'Color', cmap(1,:), 'LineWidth', linewidth);
ylim([0, 1.1*max(m1)])
ylabel('mRNA (nM)');
xlabel('Time (h)');
xlim([0, max(T)]);
xticks([0:12:max(T)]);
set(gca, 'YColor', [0 0 0], 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
yyaxis('right');
plot(T, p1, '-', 'Color', cmap(1,:), 'LineWidth', linewidth);
ylim([0, 1.1*max(p1)]);
ylabel('Protein (nM)');
xlabel('Time (h)');
xlim([0, max(T)])
xticks([0:12:max(T)]);
set(gca, 'YColor', [0 0 0], 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
legend('m_1','p_1', 'Orientation', 'horizontal', 'Location', 'southeast');

%%%%%%%% SIMULATE FIRST COMBINED CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
Y0 = zeros(25,1);
[T_0, Y_0] = odesolver(f_2_0, [- runintmax, 0], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_2_1, [- runintmax, 0], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_2_1, [0, 2*tind], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_2_2, [max(T_1) + (1e-6), max(T_1) + 6*tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

xC = Y(:, 1); mC = Y(:, 2); cC = Y(:, 3); pC = Y(:, 4); xA = Y(:, 5); mA = Y(:, 6); cA = Y(:, 7); pA = Y(:, 8); oP = Y(:, 9);
xR = Y(:,10); rR = Y(:,11); oR = Y(:,12);
x1 = Y(:,13); m1 = Y(:,14); c1 = Y(:,15); p1 = Y(:,16); x2 = Y(:,17); m2 = Y(:,18); c2 = Y(:,19); p2 = Y(:,20);
kR = Y(:,21); xF = Y(:,22); mF = Y(:,23); cF = Y(:,24); pF = Y(:,25);

figure(fplot.Number);
ax{2} = axes('Position', [0.70, 0.55, 0.2, 0.15]); box('on'); hold('on'); % pbaspect([1 1 1]);
box('on'); hold('on');
yyaxis('left');
plot(T, m1, ':', 'Color', cmap(2,:), 'LineWidth', linewidth);
ylim([0, 1.1*max(m1)])
ylabel('mRNA (nM)');
xlabel('Time (h)');
xlim([0, max(T)]);
xticks([0:12:max(T)]);
set(gca, 'YColor', [0 0 0], 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
yyaxis('right');
plot(T, p1, '-', 'Color', cmap(2,:), 'LineWidth', linewidth);
ylim([0, 1.1*max(p1)]);
ylabel('Protein (nM)');
xlabel('Time (h)');
xlim([0, max(T)])
xticks([0:12:max(T)]);
set(gca, 'YColor', [0 0 0], 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
legend('m_1','p_1', 'Orientation', 'horizontal', 'Location', 'southeast');

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
topo(1) = 1; mlim(1,:) = [0.995, 1.005]; plim(1,:) = [0.7, 1.1]; controllername{1}{1} = 'UBER-OR';
topo(2) = 2; mlim(2,:) = [0.950, 1.050]; plim(2,:) = [0.7, 1.1]; controllername{1}{2} = 'FRAG-OR';

%% %%%%% MAKE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add legend
figure(fplot.Number);
 
sp{1} = subplot(4, 4, [ 9, 10]); box('on'); hold('on'); grid('on'); pbaspect([1.5 1 1]);
sp{2} = subplot(4, 4, [11, 12]); box('on'); hold('on'); grid('on'); pbaspect([1.5 1 1]);
sp{3} = subplot(4, 4, [13, 14]); box('on'); hold('on'); pbaspect([1.5 1 1]);
sp{4} = subplot(4, 4, [15, 16]); box('on'); hold('on'); pbaspect([1.5 1 1]);
ax{3} = axes('Position', [0.35, 0.12, 0.12, 0.08]); box('on'); hold('on'); % pbaspect([1 1 1]);
ax{4} = axes('Position', [0.76, 0.12, 0.12, 0.08]); box('on'); hold('on'); % pbaspect([1 1 1]);
 
figure(fplot.Number);
for m = 1:length(topo)
    
    % make submodel
    submodel = model{topo(m)};
    
    %%%%%%%% PLOT ALL OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    %%%%%%%% SIMULATE PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    bar(ax{3}, m, Y(end, submodel.m1idx), 'FaceColor', cmap(m,:));
    set(ax{3}, 'YScale', 'log', 'YLim', [1e3 10^6.5], 'YTick', [1e3, 1e4, 1e5, 1e6], 'XTick', [], 'XLim', [0.5, 2.5]);
    ylabel(ax{3}, 'mRNA output');
    set(ax{3}, 'FontSize', 0.75*fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    plot(sp{4}, T, Y(:,submodel.p1idx)./Y_1(end,submodel.p1idx), '-', 'Color', cmap(m,:), 'LineWidth', 2);
    set(sp{4}, 'YLim', [0.7 1.1], 'XLim', [0, max(T)], 'XTick', [0:6:max(T)]);
    ylabel(sp{4}, '\Delta TL');
    xlabel(sp{4}, 'Time (h)');
    set(sp{4}, 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
    bar(ax{4}, m, Y(end, submodel.p1idx), 'FaceColor', cmap(m,:));
    set(ax{4}, 'YScale', 'log', 'YLim', [1e3 1e6], 'YTick', [1e3, 1e4, 1e5, 1e6], 'XTick', [], 'XLim', [0.5, 2.5]);
    ylabel(ax{4}, 'Protein output');
    set(ax{4}, 'FontSize', 0.75*fontsize, 'FontWeight', 'bold', 'LineWidth', outlinewidth);
    
end
