close all
clearvars


fRF = 1.57542e9;
th_in = deg2rad(90);
ph_in = deg2rad([35 62.5 90 96.5 123.5 151]);
Ps = -60; % in dBm

Pn = -85; % in dBm
LNAGain = 30;
IFGain = 50;
fIF = 4.096e6;
fLO = fRF - fIF;

fSamp = 4*fIF;          % Sample rate in Hz
Nt = 2^5;               % Number of time samples
delT = 1/fSamp;
t0 = 0;
tsamp = t0:delT:(t0+delT*(Nt-1));

Nbits = 16;
c = physconst('LightSpeed'); %propagation velocity
w = 2*pi*fRF; %center frequency
mu = 2*pi*c/w; %wavelength
r = 0.48*mu; %distance between (virtual) sensor elements
K = 40;
xPos = (-(K-1)/2:(K-1)/2).*r;
antPosIdeal = Pnt3D(xPos,0,0);

channelPhasors = ones(1,K);

arraySys = ArraySystem(antPosIdeal,channelPhasors,1,Pn,LNAGain,IFGain,fLO,fSamp,Nt,Nbits);

phaseSig = deg2rad(0);

Nsource = length(ph_in);
for ss = 1:Nsource
    S(:,ss) = PlaneWaveSignal('compExp',fRF,th_in,ph_in(ss),Ps,phaseSig);
end

portSigMat = arraySys.elements.portSignals(S,tsamp);
portSigMat_test = arraySys.elements.portSignals(S,tsamp).*lin20(LNAGain+IFGain);
[sn,si,sq] = arraySys.receiver.sigRec(portSigMat,tsamp);
% MUSIC
R = arraySys.elements.antPos.pointMatrix.';
bufferPhiD = 10;
phRange = [min(ph_in)-deg2rad(bufferPhiD),max(ph_in)+deg2rad(bufferPhiD)];
[anglesMUSIC, Z_MUSIC] = MUSICph(portSigMat,R,th_in,phRange,fRF,Nsource);

%ESPRIT
PHI = zeros(Nt,Nsource);

for k = 1:Nt

    PHI(k,:) = ESPRIT(sn(:,k),Nsource);

end
%disp(PHI)
PHI = PHI(:);
idxNaN = isnan(PHI);
PHI(idxNaN) = [];
idxInf = isinf(PHI);
PHI(idxInf) = [];

phi = log(PHI)/r;
angles = real(acos(phi.*c/(1i*w)));
angles = rad2deg(sort(angles))