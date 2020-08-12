%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all');

%% %%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

xifA = xifP;        % xifA
xirA = xirP;        % xirA
tauA = tauP;        % tauA
deltamA = deltamP;  % deltamA
betafA = betafP;    % betafA
betarA = betarP;    % betarA
gammaA = gammaP;    % gammaA
deltapA = 0;        % deltapA
alphafA = 1.62;
alpharA = 0.72;

omegagC = 500;
xifC = xifQ; xirC = xirQ; tauC = tauQ; deltamC = deltamQ;
betafC = betafQ; betarC = betarQ; gammaC = gammaQ; deltapC = deltapQ;

%% %%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tind = 120;
odeoptions = odeset();
odesolver = @ode23s;

%% %%%%% VECTORS OF BIOLOGICAL PERMISSIBLE VALUES %%%%%%%%%%%%%%%%%%%%%%%%%
v_kL = 10.^[4:7];
v_kX_const = 200; % 10.^[1:3];

v_name = { 'tetR', 'lacI', 'cI'};
v_kX =   [    350,    550,  100];
v_mu =   [    167,   0.02,   22];
v_eta =  [      1,      4,    2];

v_proteins = v_name;

v_xif = (xir1 + phiR*tau1)./v_kX;
v_xif_const = 30; % (xir1 + phiR*tau1)./v_kX_const;
v_betaf = (betar1 + gamma1)./v_kL;

% allowable plasmid copy numbers [omegagC omegagR omegagF]
allowableomega = [ ...
    [ 500, 500, 500]; ... % all genes on a high copy plasmid
    [ 100, 100, 100]; ... % all genes on med cp
    [  10,  10,  10]; ... % all genes on host
    perms([  10,  10, 500]); ...
    perms([  10, 500, 500]); ...
    perms([  10,  10, 100]); ...
    perms([  10, 100, 100]); ...
    perms([ 100, 100, 500]); ...
    perms([ 100, 500, 500]);
    ];

allowableomega = unique(allowableomega, 'rows');

%% %%%%% BIOLOIGCAL IMPLEMENTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of parameters to vary
fid = 0;
for i_1 = 1:length(allowableomega) % omegag values
    for i_2 = 1:length(v_betaf) % betafC
        for i_3 = 1:length(v_betaf) % betafA
            for i_4 = 1:length(v_betaf) % betafF
                for i_5 = 1:length(v_xif_const) % xifC
                    for i_6 = 1:length(v_xif_const) % xifA
                        for i_7 = 1:length(v_xif_const) % xifF
                            for i_8 = 1:length(v_proteins) % translational repressor
                                
                                fid = fid + 1;
                                X00(fid, 1) = allowableomega(i_1, 1);   % omegagC
                                X00(fid, 2) = allowableomega(i_1, 2);   % omegagR
                                X00(fid, 3) = allowableomega(i_1, 3);   % omegagF
                                X00(fid, 4) = v_betaf(i_2);             % betafC
                                X00(fid, 5) = v_betaf(i_3);             % betafA
                                X00(fid, 6) = v_betaf(i_4);             % betafF
                                X00(fid, 7) = v_xif_const(i_5);         % xifC
                                X00(fid, 8) = v_xif_const(i_6);         % xifA
                                X00(fid, 9) = v_xif_const(i_7);         % xifF
                                X00(fid,10) = v_xif(i_8);  % xifR
                                X00(fid,11) = v_mu(i_8);   % alpharF
                                X00(fid,12) = v_eta(i_8);  % etaF
                            end
                        end
                    end
                end
            end
        end
    end
end

% number of combinations
N = length(X00(:,1));

%% %%%%% CIRCUIT AND CONTORLLER DESIGNS TO ITTERATE OVER %%%%%%%%%%%%%%%%%%
X0(:, 1) = repmat(omegag1,[N,1]);
X0(:, 2) = xif1;
X0(:, 3) = xir1;
X0(:, 4) = tau1;
X0(:, 5) = deltam1;
X0(:, 6) = betaf1;
X0(:, 7) = betar1;
X0(:, 8) = gamma1;
X0(:, 9) = deltap1;
X0(:,10) = X00(:,1);
X0(:,11) = X00(:,7);
X0(:,12) = xirC;
X0(:,13) = tauC;
X0(:,14) = deltamC;
X0(:,15) = X00(:,4);
X0(:,16) = betarC;
X0(:,17) = gammaC;
X0(:,18) = deltapC;
X0(:,19) = X00(:,2);
X0(:,20) = X00(:,10);
X0(:,21) = xirR;
X0(:,22) = tauR;
X0(:,23) = deltarR;
X0(:,24) = rhof;
X0(:,25) = rhor;
X0(:,26) = lambda;
X0(:,27) = alphaP;
X0(:,28) = alphaR;
X0(:,29) = phiP;
X0(:,30) = phiR;
X0(:,31) = oleak;
X0(:,32) = X00(:,8);
X0(:,33) = xirA;
X0(:,34) = tauA;
X0(:,35) = deltamA;
X0(:,36) = X00(:,5);
X0(:,37) = betarA;
X0(:,38) = gammaA;
X0(:,39) = deltapA;
X0(:,40) = alphafA;
X0(:,41) = alpharA;
X0(:,42) = X00(:,3);
X0(:,43) = X00(:,9);
X0(:,44) = xirF;
X0(:,45) = tauF;
X0(:,46) = deltamF;
X0(:,47) = X00(:,6);
X0(:,48) = betarF;
X0(:,49) = gammaF;
X0(:,50) = deltapF;
X0(:,51) = alphafF;
X0(:,52) = X00(:,11);
X0(:,53) = X00(:,12);

%% %%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('simulationsettings');

parfor i = 1:length(X0(:,1))
    
    disp(['... ',num2str(i)]); warning('off');
    
    % get i th parameter set
    X = X0(i,:);
    
    f_1 = @(T,Y) frag_or(T, Y, [X(2:9), X(2:9)], [X(11:18)], [X(20:25)], [X(26:31)], [    0,    0, X(10), X(19), X(42)], [X(32:41)], [X(43:53)]);
    f_2 = @(T,Y) frag_or(T, Y, [X(2:9), X(2:9)], [X(11:18)], [X(20:25)], [X(26:31)], [ X(1),    0, X(10), X(19), X(42)], [X(32:41)], [X(43:53)]);
    f_3 = @(T,Y) frag_or(T, Y, [X(2:9), X(2:9)], [X(11:18)], [X(20:25)], [X(26:31)], [ X(1), X(1), X(10), X(19), X(42)], [X(32:41)], [X(43:53)]);
    
    % simulation
    Y0 = zeros(25,1);
    [T_1, Y_1] = odesolver(f_1, [0, tind], Y0, odeoptions);
    Y0 = Y_1(end,:);
    [T_2, Y_2] = odesolver(f_2, [tind + (1e-6), 2*tind], Y0, odeoptions);
    Y0 = Y_2(end,:);
    [T_3, Y_3] = odesolver(f_3, [2*tind + (1e-6), 3*tind], Y0, odeoptions);
    
    % save results to a text file with high precision
    data_0 = [T_1, Y_1];
    data_1 = [T_2, Y_2];
    data_2 = [T_3, Y_3];
    
    dlmwrite(['i_',num2str(i),'_dynamics_0.txt'], data_0, 'precision', '%.9f');
    dlmwrite(['i_',num2str(i),'_dynamics_1.txt'], data_1, 'precision', '%.9f');
    dlmwrite(['i_',num2str(i),'_dynamics_2.txt'], data_2, 'precision', '%.9f');
    
end

