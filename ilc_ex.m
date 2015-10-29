% Input parameters
fn = 1;
zeta = 0.1;
kp = 0;
snr = inf;
cutOffFrq = inf;
resetIc = 1;

fd = 1000;
numOfPer = 5;

numOfIter = 100000;

% System
wn = fn*2*pi;
s = tf('s');
sys = wn^2/(s^2 + 2*zeta*wn*s + wn^2);
sys_ss = ss(sys);


% Desired Trajectory
Td = 1/fd;
yd =@(tt) (mod(tt,Td) <= Td/2).*((4/Td)*mod(tt,Td) - 1) ...
        + (mod(tt,Td) > Td/2).*(3 - (4/Td)*mod(tt,Td));

% Time vecotor
fs = 1000*fd;
tt = 0:(1/fs):(2*numOfPer*(1/fd)-(1/fs));
tt = tt(:);

L = numel(tt);
LL = L/2;

f = (-LL/2:LL/2 - 1)*fs/LL;
f = f(:);

yd = yd(tt);

%[uff,Ginv] = miifc_lti(sys_ss,tt,yd,snr,cutOffFrq,resetIc,numOfIter);
[Emax,Erms] = ifc_lti(sys_ss,tt,yd,snr,cutOffFrq,resetIc,numOfIter);


figure; hold all;
  plot(Emax,'--k');
  plot(Erms,'r');
  set(gca,'xscale','log')
