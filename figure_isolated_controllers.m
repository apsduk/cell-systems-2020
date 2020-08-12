%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all');

tind = 36;
tstart = 12;
tmax = 120;

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

% ornap dissociation constant
kX0 = 200;

% calculate xif1 and xif2
xif1 = (xir1 + phiR*tau1)./kX0;
xif2 = xif1;

% make circuit vectors
% betaf2 = 0;
circPR = [xif1 xir1 tau1 deltam1 betaf1 betar1 gamma1 deltap1 xif2 xir2 tau2 deltam2 betaf2 betar2 gamma2 deltap2];

%% %%%%% TRANSCRIPTIONAL CONTROLLER 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specific circuit parameters
oleak = 1;
omegagP = 100;
omegagQ = 10;
betafP = 100.03;
betafQ = 100.03;
xifP = 55; % o-rnap
xifQ = 55; % o-rnap
alpharQ = 167;
etaQ = 1;

% make vectors
ornapPR = [xifP xirP tauP deltamP betafP betarP gammaP deltaoP];
rnapcontPR = [xifQ xirQ tauQ deltamQ betafQ betarQ gammaQ deltapQ alphafQ alpharQ etaQ];
hostPR = [lambda alphaP alphaR phiP phiR oleak];

% ode model
f_0 = @(T,Y) uber(T, Y, circPR, ornapPR, [], hostPR, [omegagP       0       0 omegagQ], rnapcontPR, []);
f_1 = @(T,Y) uber(T, Y, circPR, ornapPR, [], hostPR, [omegagP omegag1       0 omegagQ], rnapcontPR, []);
f_2 = @(T,Y) uber(T, Y, circPR, ornapPR, [], hostPR, [omegagP omegag1 omegag2 omegagQ], rnapcontPR, []);

% simulation
Y0 = zeros(17,1);
[T_0, Y_0] = odesolver(f_0, [0, tmax], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_1, [0, tmax], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_1, [0, tstart], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_2, [tstart + (1e-6), tstart + tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

% species dynamics
xP = Y(:, 1); mP = Y(:, 2); cP = Y(:, 3); oP = Y(:, 4);
x1 = Y(:, 5); m1 = Y(:, 6); c1 = Y(:, 7); p1 = Y(:, 8);
x2 = Y(:, 9); m2 = Y(:,10); c2 = Y(:,11); p2 = Y(:,12);
kP = Y(:,13); xQ = Y(:,14); mQ = Y(:,15); cQ = Y(:,16); pQ = Y(:,17);

% save results
isolatedcontroller{1}.T = T;
isolatedcontroller{1}.m1 = m1;
isolatedcontroller{1}.m2 = m2;
isolatedcontroller{1}.p1 = p1;
isolatedcontroller{1}.p2 = p2;

%% %%%%% TRANSCRIPTIONAL CONTROLLER 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xifA = xifP; xirA = xirP; tauA = tauP; deltamA = deltamP; betafA = betafP; betarA = betarP; gammaA = gammaP; deltapA = 0;
xifC = xifP; xirC = xirP; tauC = tauP; deltamC = deltamP; betafC = betafP; betarC = betarP; gammaC = gammaP;

% specific circuit parameters
oleak = 1;
xifC = 30; % h-rnap
xifA = 30; % h-rnap
betafC = 100.03;
betafA = 100.03;
alphafA = 1.62;  %  :: From Segall-Shapiro2014 ka = 4.5e5 M/s
alpharA = 0.72;  %  :: From Segall-Shapiro2014 kd = 2.0e-4 1/s
omegagC = 500;

% make vectors
ornapPR = [xifC xirC tauC deltamC betafC betarC gammaC deltaoP];
rnapcontPR = [xifA xirA tauA deltamA betafA betarA gammaA deltapA alphafA alpharA];
hostPR = [lambda alphaP alphaR phiP phiR oleak];

% ode model
f_0 = @(T,Y) frag(T, Y, circPR, ornapPR, [], hostPR, [      0       0 omegagC], rnapcontPR, []);
f_1 = @(T,Y) frag(T, Y, circPR, ornapPR, [], hostPR, [omegag1       0 omegagC], rnapcontPR, []);
f_2 = @(T,Y) frag(T, Y, circPR, ornapPR, [], hostPR, [omegag1 omegag2 omegagC], rnapcontPR, []);

% simulation
Y0 = zeros(17,1);
[T_0, Y_0] = odesolver(f_0, [0, tmax], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_1, [0, tmax], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_1, [0, tstart], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_2, [tstart + (1e-6), tstart + tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

% species dynamics
xC = Y(:, 1); mC = Y(:, 2); cC = Y(:, 3); pC = Y(:, 4);
xA = Y(:, 5); mA = Y(:, 6); cA = Y(:, 7); pA = Y(:, 8); oP = Y(:, 9);
x1 = Y(:,10); m1 = Y(:,11); c1 = Y(:,12); p1 = Y(:,13);
x2 = Y(:,14); m2 = Y(:,15); c2 = Y(:,16); p2 = Y(:,17);

% save results
isolatedcontroller{2}.T = T;
isolatedcontroller{2}.m1 = m1;
isolatedcontroller{2}.m2 = m2;
isolatedcontroller{2}.p1 = p1;
isolatedcontroller{2}.p2 = p2;

%% %%%%% TRANSLATIONAL CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specific circuit parameters
xif1 = 30; % h-rnap
xif2 = 30; % h-rnap

% specific controller parameters
omegagR = 100;
omegagF = 10;
xifR = 20;
alpharF = 0.02;
etaF = 4;
xifF = 30; % h-rnap
betafF = 100.03;

% make circuit vectors
circPR = [xif1 xir1 tau1 deltam1 betaf1 betar1 gamma1 deltap1 xif2 xir2 tau2 deltam2 betaf2 betar2 gamma2 deltap2];

% make vectors
oriboPR = [xifR xirR tauR deltarR rhof rhor];
ribocontPR = [xifF xirF tauF deltamF betafF betarF gammaF deltapF alphafF alpharF etaF];
hostPR = [lambda alphaP alphaR phiP phiR];

% ode model
f_0 = @(T,Y) translationalcontroller(T, Y, circPR, [], oriboPR, hostPR, [omegagR       0       0 omegagF], [], ribocontPR);
f_1 = @(T,Y) translationalcontroller(T, Y, circPR, [], oriboPR, hostPR, [omegagR omegag1       0 omegagF], [], ribocontPR);
f_2 = @(T,Y) translationalcontroller(T, Y, circPR, [], oriboPR, hostPR, [omegagR omegag1 omegag2 omegagF], [], ribocontPR);

% simulation
Y0 = zeros(16,1);
[T_0, Y_0] = odesolver(f_0, [0, tmax], Y0, odeoptions);
[T_0, Y_0] = odesolver(f_1, [0, tmax], Y_0(end,:), odeoptions);
[T_1, Y_1] = odesolver(f_1, [0, tstart], Y_0(end,:), odeoptions);
[T_2, Y_2] = odesolver(f_2, [tstart + (1e-6), tstart + tind], Y_1(end,:), odeoptions);
T = [T_1; T_2];
Y = [Y_1; Y_2];

% species dynamics
xR = Y(:, 1); rR = Y(:, 2); oR = Y(:, 3);
x1 = Y(:, 4); m1 = Y(:, 5); c1 = Y(:, 6); p1 = Y(:, 7);
x2 = Y(:, 8); m2 = Y(:, 9); c2 = Y(:,10); p2 = Y(:,11);
kR = Y(:,12); xF = Y(:,13); mF = Y(:,14); cF = Y(:,15); pF = Y(:,16);

% save results
isolatedcontroller{3}.T = T;
isolatedcontroller{3}.m1 = m1;
isolatedcontroller{3}.m2 = m2;
isolatedcontroller{3}.p1 = p1;
isolatedcontroller{3}.p2 = p2;

%% %%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [1, 2; 3, 4; 5, 6];

linewidth = 1;

fplot = figure;
fplot.Units = 'centimeters';
fplot.Position = [5, 5, 15, 20];

figure(fplot.Number);
for i = 1:length(isolatedcontroller)
    subplot(3, 2, fid(i,1)); box('on'); hold('on');
    plot(isolatedcontroller{i}.T, isolatedcontroller{i}.m1, '-', 'LineWidth', linewidth);
    plot(isolatedcontroller{i}.T, isolatedcontroller{i}.m2, '-', 'LineWidth', linewidth);
    ylabel('mRNA (nM)');
    xlabel('Time (h)');
    xlim([0, max(T)]);
    ylim([0, 1.1*max(isolatedcontroller{i}.m1)])
    xticks([0:12:max(T)]);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    if betaf2 == 0
        legend('m_1','sRNA','Location','southeast');
    else
        legend('m_1','m_2','Location','southeast'); % ,'Orientation','horizontal');
    end
    pbaspect([1, 1, 1]);
    subplot(3, 2, fid(i,2)); box('on'); hold('on');
    plot(isolatedcontroller{i}.T, isolatedcontroller{i}.p1, '-', 'LineWidth', linewidth);
    if betaf2 ~= 0
        plot(isolatedcontroller{i}.T, isolatedcontroller{i}.p2, '-', 'LineWidth', linewidth);
    end
    xlabel('Time (h)');
    xlim([0, max(T)])
    xticks([0:12:max(T)]);
    ylim([0, 1.1*max(isolatedcontroller{i}.p1)]);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    if betaf2 == 0
        legend('p_1','Location','southeast');
    else
        legend('p_1','p_2','Location','southeast'); % ,'Orientation','horizontal');
    end
    ylabel('Protein (nM)');
    pbaspect([1, 1, 1]);
    
end
