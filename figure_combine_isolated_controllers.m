%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all');

tind = 36;
tstart = 12;
tmax = 300;

odesolver = @ode23s;
odeoptions = odeset();

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

kX0 = 200;

xif1 = (xir1 + phiR*tau1)./kX0;
xif2 = xif1;

% make circuit vectors
% betaf2 = 0;
circPR = [xif1 xir1 tau1 deltam1 betaf1 betar1 gamma1 deltap1 xif2 xir2 tau2 deltam2 betaf2 betar2 gamma2 deltap2];

% host prameters with o leak
hostPR = [lambda alphaP alphaR phiP phiR oleak];

%% %%%%% TRANSCRIPTIONAL CONTROLLER 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% %%%%% TRANSCRIPTIONAL CONTROLLER 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%% TRANSLATIONAL CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% %%%%% MAKE ODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_1_0 = @(T,Y) uber_or(T, Y, circPR, ornapPR1, oriboPR, hostPR, [omegagP omegagR       0       0 omegagF omegagQ], rnapcontPR1, ribocontPR);
f_1_1 = @(T,Y) uber_or(T, Y, circPR, ornapPR1, oriboPR, hostPR, [omegagP omegagR omegag1       0 omegagF omegagQ], rnapcontPR1, ribocontPR);
f_1_2 = @(T,Y) uber_or(T, Y, circPR, ornapPR1, oriboPR, hostPR, [omegagP omegagR omegag1 omegag2 omegagF omegagQ], rnapcontPR1, ribocontPR);

f_2_0 = @(T,Y) frag_or(T, Y, circPR, ornapPR2, oriboPR, hostPR, [      0       0 omegagC omegagR omegagF], rnapcontPR2, ribocontPR);
f_2_1 = @(T,Y) frag_or(T, Y, circPR, ornapPR2, oriboPR, hostPR, [omegag1       0 omegagC omegagR omegagF], rnapcontPR2, ribocontPR);
f_2_2 = @(T,Y) frag_or(T, Y, circPR, ornapPR2, oriboPR, hostPR, [omegag1 omegag2 omegagC omegagR omegagF], rnapcontPR2, ribocontPR);

%% %%%%% SIMULATE UBER-OR COMBINED CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
Y0 = zeros(25,1);
[T_0, Y_0] = odesolver(f_1_0, [0, tmax], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_1_1, [0, tmax], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_1_1, [0, tstart], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_1_2, [tstart + (1e-6), tstart + tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

xP = Y(:, 1); mP = Y(:, 2); cP = Y(:, 3); oP = Y(:, 4);
xR = Y(:, 5); rR = Y(:, 6); oR = Y(:, 7);
x1 = Y(:, 8); m1 = Y(:, 9); c1 = Y(:,10); p1 = Y(:,11);
x2 = Y(:,12); m2 = Y(:,13); c2 = Y(:,14); p2 = Y(:,15);
kR = Y(:,16); xF = Y(:,17); mF = Y(:,18); cF = Y(:,19); pF = Y(:,20);
kP = Y(:,21); xQ = Y(:,22); mQ = Y(:,23); cQ = Y(:,24); pQ = Y(:,25);

% save results
dualcontroller{1}.T = T;
dualcontroller{1}.m1 = m1;
dualcontroller{1}.m2 = m2;
dualcontroller{1}.p1 = p1;
dualcontroller{1}.p2 = p2;

%% %%%%% SIMULATE FRAG-OR COMBINED CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
Y0 = zeros(25,1);
[T_0, Y_0] = odesolver(f_2_0, [0, tmax], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_2_1, [0, tmax], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_2_1, [0, tstart], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_2_2, [tstart + (1e-6), tstart + tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

xC = Y(:, 1); mC = Y(:, 2); cC = Y(:, 3); pC = Y(:, 4); xA = Y(:, 5); mA = Y(:, 6); cA = Y(:, 7); pA = Y(:, 8); oP = Y(:, 9);
xR = Y(:,10); rR = Y(:,11); oR = Y(:,12);
x1 = Y(:,13); m1 = Y(:,14); c1 = Y(:,15); p1 = Y(:,16); x2 = Y(:,17); m2 = Y(:,18); c2 = Y(:,19); p2 = Y(:,20);
kR = Y(:,21); xF = Y(:,22); mF = Y(:,23); cF = Y(:,24); pF = Y(:,25);

% save results
dualcontroller{2}.T = T;
dualcontroller{2}.m1 = m1;
dualcontroller{2}.m2 = m2;
dualcontroller{2}.p1 = p1;
dualcontroller{2}.p2 = p2;

%% %%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [1, 2; 3, 4];

linewidth = 1.5;

fplot = figure;
fplot.Units = 'centimeters';
fplot.Position = [5, 5, 15, 20];

figure(fplot.Number);
for i = 1:length(dualcontroller)
    subplot(2, 2, fid(i,1)); box('on'); hold('on');
    plot(dualcontroller{i}.T, dualcontroller{i}.m1, '-', 'LineWidth', linewidth);
    plot(dualcontroller{i}.T, dualcontroller{i}.m2, '-', 'LineWidth', linewidth);
    ylabel('mRNA (nM)');
    xlabel('Time (h)');
    xlim([0, max(T)]);
    ylim([0, 1.1*max(dualcontroller{i}.m1)])
    xticks([0:12:max(T)]);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    legend('m_1','m_2','Location','southeast'); % ,'Orientation','horizontal');
    subplot(2, 2, fid(i,2)); box('on'); hold('on');
    plot(dualcontroller{i}.T, dualcontroller{i}.p1, '-', 'LineWidth', linewidth);
    plot(dualcontroller{i}.T, dualcontroller{i}.p2, '-', 'LineWidth', linewidth);
    xlabel('Time (h)');
    xlim([0, max(T)])
    xticks([0:12:max(T)]);
    ylim([0, 1.1*max(dualcontroller{i}.p1)]);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    legend('p_1','p_2','Location','southeast'); % ,'Orientation','horizontal');
    ylabel('Protein (nM)');

end