function dYdt = frag(T, Y, circPR, ornapPR, ~, hostPR, plasmidPR, rnapcontPR, ~)

%% %%%%% SPECIES CURRENT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xC = Y( 1); mC = Y( 2); cC = Y( 3); pC = Y( 4);
xA = Y( 5); mA = Y( 6); cA = Y( 7); pA = Y( 8); oP = Y( 9);
x1 = Y(10); m1 = Y(11); c1 = Y(12); p1 = Y(13);
x2 = Y(14); m2 = Y(15); c2 = Y(16); p2 = Y(17);

%% %%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% circuit parameters
xif1 = circPR( 1); xir1 = circPR( 2); tau1 = circPR( 3); deltam1 = circPR( 4); betaf1 = circPR( 5); betar1 = circPR( 6); gamma1 = circPR( 7); deltap1 = circPR( 8);
xif2 = circPR( 9); xir2 = circPR(10); tau2 = circPR(11); deltam2 = circPR(12); betaf2 = circPR(13); betar2 = circPR(14); gamma2 = circPR(15); deltap2 = circPR(16);

% core o-rnap
xifC = ornapPR(1); xirC = ornapPR(2); tauC = ornapPR(3); deltamC = ornapPR(4); betafC = ornapPR(5); betarC = ornapPR(6); gammaC = ornapPR(7); deltapC = ornapPR(8);

% alpha subunit o-rnap controller
xifA = rnapcontPR(1); xirA = rnapcontPR(2); tauA = rnapcontPR(3); deltamA = rnapcontPR(4); betafA = rnapcontPR(5); betarA = rnapcontPR(6); gammaA = rnapcontPR(7); deltapA = rnapcontPR(8); alphafA = rnapcontPR(9); alpharA = rnapcontPR(10);

% host parameters
lambda = hostPR(1); alphaP = hostPR(2); alphaR = hostPR(3); phiP = hostPR(4); phiR = hostPR(5); oleak = hostPR(6);

% plasmid parameters
omegag1 = plasmidPR(1); omegag2 = plasmidPR(2); omegagC = plasmidPR(3); omegagA = omegag1 + omegag2;

%% %%%%% ACCOUNT FOR MULTIPLE RESOURCE USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omegag1 = phiP*omegag1; omegag2 = phiP*omegag2;
tau1 = phiR*tau1; tau2 = phiR*tau2;
omegagC = phiP*omegagC; omegagA = phiP*omegagA;
tauC = phiR*tauC; tauA = phiR*tauA;

%% %%%%% APPLY CONSERVATION LAWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gA = omegagA/lambda - xA;
gC = omegagC/lambda - xC;
g1 = omegag1/lambda - x1;
g2 = omegag2/lambda - x2;

Ptot = alphaP/lambda;
Rtot = alphaR/lambda;

hP = Ptot - xC - xA;
hR = Rtot - cA - cC - c1 - c2;

%% %%%%% ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% core o-RNAP
dxC = xifC*gC*hP - (xirC + tauC)*xC - lambda*xC;
dmC = tauC*xC - betafC*mC*hR + (betarC + gammaC)*cC - (deltamC + lambda)*mC;
dcC = betafC*mC*hR - (betarC + gammaC)*cC - lambda*cC;
dpC = gammaC*cC - (deltapC + lambda)*pC ...
    - alphafA*pA*pC + alpharA*oP;

% alpha o-RNAP
dxA = xifA*gA*hP - (xirA + tauA)*xA - lambda*xA;
dmA = tauA*xA - betafA*mA*hR + (betarA + gammaA)*cA - (deltamA + lambda)*mA;
dcA = betafA*mA*hR - (betarA + gammaA)*cA - lambda*cA;
dpA = gammaA*cA - (deltapA + lambda)*pA ...
    - alphafA*pA*pC + alpharA*oP;

% orthohonal RNAP
doP = oleak + alphafA*pA*pC - alpharA*oP ...
    - (deltapC + lambda)*oP ...
    - xif1*g1*oP + (xir1 + tau1)*x1 ...
    - xif2*g2*oP + (xir2 + tau2)*x2;

% circuit ODEs
dx1 = xif1*g1*oP - (xir1 + tau1)*x1 - lambda*x1;
dm1 = tau1*x1 - betaf1*m1*hR + (betar1 + gamma1)*c1 - (deltam1 + lambda)*m1;
dc1 = betaf1*m1*hR - (betar1 + gamma1)*c1 - lambda*c1;
dp1 = gamma1*c1 - (deltap1 + lambda)*p1;
dx2 = xif2*g2*oP - (xir2 + tau2)*x2 - lambda*x2;
dm2 = tau2*x2 - betaf2*m2*hR + (betar2 + gamma2)*c2 - (deltam2 + lambda)*m2;
dc2 = betaf2*m2*hR - (betar2 + gamma2)*c2 - lambda*c2;
dp2 = gamma2*c2 - (deltap2 + lambda)*p2;

%% %%%%% RETURN dYdt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dYdt = [dxC; dmC; dcC; dpC; dxA; dmA; dcA; dpA; doP; dx1; dm1; dc1; dp1; dx2; dm2; dc2; dp2];

end