function dYdt = translationalcontroller(T, Y, circPR, ~, oriboPR, hostPR, plasmidPR, ~, ribocontPR)

%% %%%%% SPECIES CURRENT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xR = Y( 1); rR = Y( 2); oR = Y( 3);
x1 = Y( 4); m1 = Y( 5); c1 = Y( 6); p1 = Y( 7);
x2 = Y( 8); m2 = Y( 9); c2 = Y(10); p2 = Y(11);
kR = Y(12); xF = Y(13); mF = Y(14); cF = Y(15); pF = Y(16);

%% %%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% circuit parameters
xif1 = circPR( 1); xir1 = circPR( 2); tau1 = circPR( 3); deltam1 = circPR( 4); betaf1 = circPR( 5); betar1 = circPR( 6); gamma1 = circPR( 7); deltap1 = circPR( 8);
xif2 = circPR( 9); xir2 = circPR(10); tau2 = circPR(11); deltam2 = circPR(12); betaf2 = circPR(13); betar2 = circPR(14); gamma2 = circPR(15); deltap2 = circPR(16);

% o-ribosome parameters
xifR = oriboPR(1); xirR = oriboPR(2); tauR = oriboPR(3); deltarR = oriboPR(4); rhof = oriboPR(5); rhor = oriboPR(6);

% o-ribosome controller
xifF = ribocontPR(1); xirF = ribocontPR(2); tauF = ribocontPR(3); deltamF = ribocontPR(4); betafF = ribocontPR(5); betarF = ribocontPR(6); gammaF = ribocontPR(7); deltapF = ribocontPR(8); alphafF = ribocontPR(9); alpharF = ribocontPR(10); etaF = ribocontPR(11);

% host parameters
lambda = hostPR(1); alphaP = hostPR(2); alphaR = hostPR(3); phiP = hostPR(4); phiR = hostPR(5);

% plasmid copy numbers
omegagR = plasmidPR(1); omegag1 = plasmidPR(2); omegag2 = plasmidPR(3); omegagF = plasmidPR(4);

%% %%%%% ACCOUNT FOR MULTIPLE RESOURCE USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omegagR = phiP*omegagR; omegag1 = phiP*omegag1; omegag2 = phiP*omegag2; omegagF = phiP*omegagF;
tau1 = phiR*tau1; tau2 = phiR*tau2; tauF = phiR*tauF;

%% %%%%% APPLY CONSERVATION LAWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gR = omegagR - xR - kR;
g1 = omegag1 - x1;
g2 = omegag2 - x2;
gF = omegagF - xF;

Ptot = alphaP/lambda;
Rtot = alphaR/lambda;

hP = Ptot - xR - xF - x1 - x2;
hR = Rtot - oR - c1 - c2 - cF;

%% %%%%% ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orthogonal ribosomes
dxR = xifR*gR*hP - (xirR + tauR)*xR - lambda*xR;
drR = tauR*xR - (deltarR + lambda)*rR - rhof*rR*hR + rhor*oR;
doR = rhof*rR*hR - rhor*oR - lambda*oR ...
    - betaf1*m1*oR + (betar1 + gamma1)*c1 ...
    - betaf2*m2*oR + (betar2 + gamma2)*c2 ...
    - betafF*mF*oR + (betarF + gammaF)*cF;

% circuit ODEs
dx1 = xif1*g1*hP - (xir1 + tau1)*x1 - lambda*x1;
dm1 = tau1*x1 - betaf1*m1*oR + (betar1 + gamma1)*c1 - (deltam1 + lambda)*m1;
dc1 = betaf1*m1*oR - (betar1 + gamma1)*c1 - lambda*c1;
dp1 = gamma1*c1 - (deltap1 + lambda)*p1;
dx2 = xif2*g2*hP - (xir2 + tau2)*x2 - lambda*x1;
dm2 = tau2*x2 - betaf2*m2*oR + (betar2 + gamma2)*c2 - (deltam2 + lambda)*m2;
dc2 = betaf2*m2*oR - (betar2 + gamma2)*c2 - lambda*c2;
dp2 = gamma2*c2 - (deltap2 + lambda)*p2;

% o-ribosome controller protein
dkR = alphafF*gR*pF^etaF - alpharF*kR - lambda*kR;
dxF = xifF*gF*hP - (xirF + tauF)*xF - lambda*xF;
dmF = tauF*xF - betafF*mF*oR + (betarF + gammaF)*cF - (deltamF + lambda)*mF;
dcF = betafF*mF*oR - (betarF + gammaF)*cF - lambda*cF;
dpF = gammaF*cF - (deltapF + lambda)*pF - etaF*alphafF*gR*pF^etaF + etaF*alpharF*kR;

%% %%%%% RETURN dYdt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dYdt = [dxR; drR; doR; dx1; dm1; dc1; dp1; dx2; dm2; dc2; dp2; dkR; dxF; dmF; dcF; dpF];

end