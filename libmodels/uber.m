function dYdt = uber(T, Y, circPR, ornapPR, ~, hostPR, plasmidPR, rnapcontPR, ~)

%% %%%%% SPECIES CURRENT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xP = Y( 1); mP = Y( 2); cP = Y( 3); oP = Y( 4);
x1 = Y( 5); m1 = Y( 6); c1 = Y( 7); p1 = Y( 8);
x2 = Y( 9); m2 = Y(10); c2 = Y(11); p2 = Y(12);
kP = Y(13); xQ = Y(14); mQ = Y(15); cQ = Y(16); pQ = Y(17);

%% %%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% circuit parameters
xif1 = circPR( 1); xir1 = circPR( 2); tau1 = circPR( 3); deltam1 = circPR( 4); betaf1 = circPR( 5); betar1 = circPR( 6); gamma1 = circPR( 7); deltap1 = circPR( 8);
xif2 = circPR( 9); xir2 = circPR(10); tau2 = circPR(11); deltam2 = circPR(12); betaf2 = circPR(13); betar2 = circPR(14); gamma2 = circPR(15); deltap2 = circPR(16);

% o-rnap
xifP = ornapPR(1); xirP = ornapPR(2); tauP = ornapPR(3); deltamP = ornapPR(4); betafP = ornapPR(5); betarP = ornapPR(6); gammaP = ornapPR(7); deltaoP = ornapPR(8);

% o-rnap controller
xifQ = rnapcontPR(1); xirQ = rnapcontPR(2); tauQ = rnapcontPR(3); deltamQ = rnapcontPR(4); betafQ = rnapcontPR(5); betarQ = rnapcontPR(6); gammaQ = rnapcontPR(7); deltapQ = rnapcontPR(8); alphafQ = rnapcontPR(9); alpharQ = rnapcontPR(10); etaQ = rnapcontPR(11);

% host parameters
lambda = hostPR(1); alphaP = hostPR(2); alphaR = hostPR(3); phiP = hostPR(4); phiR = hostPR(5); oleak = hostPR(6);

% plasmid copy numbers production
omegagP = plasmidPR(1); omegag1 = plasmidPR(2); omegag2 = plasmidPR(3); omegagQ = plasmidPR(4);

%% %%%%% ACCOUNT FOR MULTIPLE RESOURCE USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omegagP = phiP*omegagP; omegag1 = phiP*omegag1; omegag2 = phiP*omegag2; omegagQ = phiP*omegagQ;
tau1 = phiR*tau1; tau2 = phiR*tau2; tauP = phiR*tauP; tauQ = phiR*tauQ;

%% %%%%% APPLY CONSERVATION LAWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gP = omegagP/lambda - xP - kP;
g1 = omegag1/lambda - x1;
g2 = omegag2/lambda - x2;
gQ = omegagQ/lambda - xQ;

Ptot = alphaP/lambda;
Rtot = alphaR/lambda;

hP = Ptot;
hR = Rtot - cP - c1 - c2 - cQ;

%% %%%%% ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orthohonal RNAP
dxP = + xifP*gP*oP - (xirP + tauP)*xP - lambda*xP;
dmP = tauP*xP - betafP*mP*hR + (betarP + gammaP)*cP - (deltamP + lambda)*mP;
dcP = betafP*mP*hR - (betarP + gammaP)*cP - lambda*cP;
doP = oleak + gammaP*cP - (deltaoP + lambda)*oP...
    - xifP*gP*oP + (xirP + tauP)*xP ...
    - xif1*g1*oP + (xir1 + tau1)*x1 ...
    - xif2*g2*oP + (xir2 + tau2)*x2 ...
    - xifQ*gQ*oP + (xirQ + tauQ)*xQ;

% circuit ODEs
dx1 = xif1*g1*oP - (xir1 + tau1)*x1 - lambda*x1;
dm1 = tau1*x1 - betaf1*m1*hR + (betar1 + gamma1)*c1 - (deltam1 + lambda)*m1;
dc1 = betaf1*m1*hR - (betar1 + gamma1)*c1 - lambda*c1;
dp1 = gamma1*c1 - (deltap1 + lambda)*p1;
dx2 = xif2*g2*oP - (xir2 + tau2)*x2 - lambda*x1;
dm2 = tau2*x2 - betaf2*m2*hR + (betar2 + gamma2)*c2 - (deltam2 + lambda)*m2;
dc2 = betaf2*m2*hR - (betar2 + gamma2)*c2 - lambda*c2;
dp2 = gamma2*c2 - (deltap2 + lambda)*p2;

% o-rnap controller
dkP = alphafQ*gP*pQ^etaQ - alpharQ*kP - lambda*kP;
dxQ = + xifQ*gQ*oP - (xirQ + tauQ)*xQ - lambda*xQ;
dmQ = tauQ*xQ - betafQ*mQ*hR + (betarQ + gammaQ)*cQ - (deltamQ + lambda)*mQ;
dcQ = betafQ*mQ*hR - (betarQ + gammaQ)*cQ - lambda*cQ;
dpQ = gammaQ*cQ - (deltapQ + lambda)*pQ - etaQ*alphafQ*gP*pQ^etaQ + etaQ*alpharQ*kP;

%% %%%%% RETURN dYdt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dYdt = [dxP; dmP; dcP; doP; dx1; dm1; dc1; dp1; dx2; dm2; dc2; dp2; dkP; dxQ; dmQ; dcQ; dpQ];

end