function [uff, Ginv] = miifc_lti(sys, tt, yd, snr,  ...
                                cutOffFreq, resestIc, numOfIter)
% Modeling-free inversion-based iterative feedforward control (miifc)
%
% Algorithm from Kim and Zou 2013 (IEEE/ASME Trans. on Mech.)
%
%function [uff, Ginv] = miifc_lti(sys, tt, yd, snr,  ...
%                                cutOffFreq, resestIc, numOfIter)
%

fs = 1/(tt(2) - tt(1));

% Column vectors
tt = tt(:);
yd = yd(:);

% Lenght of time vecotor (and half length);
L = numel(tt);
LL = L/2;

% Error Metrics
Emax =@(yd, y) norm(yd - y,inf)./norm(yd,inf) * 100;
Erms =@(yd, y) norm(yd - y,2)./norm(yd,2) * 100;

% Frequency vector
f = (-LL/2:LL/2 - 1)*fs/LL;
f = f(:);

% Start index for fft (half of the yd)
si = LL + 1;

% FFT fo desired traj.

Yd = fftshift(fft(yd(si:end)));
Yd = Yd.*(abs(f) < cutOffFreq);

% Find peaks
[val,index] = findpeaks(abs(Yd),'SORTSTR','descend');
index(val < 0.01) = [];

% Truncated
yd_truc = ifft(ifftshift(Yd));
yd_truc = yd_truc(:);
yd_truc = repmat(yd_truc,2,1);

% Desired Trajectory
figure; hold all;
  title('Desired Trajectory','fontsize',20);
  plot(tt,yd,'k')
  plot(tt,yd_truc,'--r');
  xlabel('Time (s)','fontsize',20);
  grid on;

% Bode
frange = logspace(-4,4,1000);
[mag, phase] = bode(sys,frange*2*pi);
mag = mag(:);
phase = phase(:);

fq = abs(f(index(1:2:end)));
[mag_yd, phase_yd] = bode(sys,fq*2*pi);
mag_yd = mag_yd(:);
phase_yd = phase_yd(:);

figure;
  subplot(211); hold all;
    semilogx(frange,mag2db(mag),'k');
    semilogx(fq,mag2db(mag_yd),'or');
    set(gca,'xscale','log');
    ylabel('Mag (dB)','fontsize',20);
    grid on
  subplot(212); hold all;
    semilogx(frange,phase,'k');
    semilogx(fq,phase_yd,'or');
    set(gca,'xscale','log');
    ylabel('Phase (deg)','fontsize',20);
    xlabel('Frequency (Hz)','fontsize',20);
    grid on

% ---- Simulation ---- %

Uff = 0.00001.*Yd;
uff = ifft(ifftshift(Uff));
uff = repmat(uff,2,1);

Emax_vector = zeros(1,numOfIter);
Erms_vector = zeros(1,numOfIter);

figure;
   subplot(511), hold all;
    h_title = title('Iteration = ','fontsize',20);
    plot(tt,yd_truc,'r');
    h_y = plot(tt,zeros(size(tt)));
    h_err = plot(tt,zeros(size(tt)));
    ylim([-2,2]);
    grid on

  subplot(512); hold all;
    h_uff = plot(tt,uff,'m');
    xlabel('Time (s)','fontsize',20);
    grid on

  subplot(513); hold all;
    h_Ginv = plot(f,0*f);
    xlim([0,max(f)])
    xlim([0,cutOffFreq + 1])
    grid on

  subplot(514); hold all;
    h_Uff = plot(f,0*f);
    xlim([0,cutOffFreq + 1])
    xlabel('Frequency (Hz)','fontsize',20);
    grid on

  subplot(515); hold all
     h_er_max = plot(0,0,'k');
     h_er_rms = plot(0,0,'r');
     xlabel('Iteration (n)','fontsize',20);
     grid on

  set(gcf,'Position',[20,20,1000,1000]);

x0 = [0; 0];
for i=1:numOfIter
  tic
  [y,~,x] = lsim(sys,uff,tt,x0);
  y = awgn(y,snr,'measured');
  toc

  if resestIc
    x0 = [0; 0];
  else
    x0 = x(end,:);
  end

  Erms_vector(i) = Erms(yd_truc(si:end),y(si:end));
  Emax_vector(i) = Emax(yd_truc(si:end),y(si:end));

  % Update
  Yk = fftshift(fft(y(si:end)));

  Ginv = (Uff./Yk);
  Ginv(isnan(Ginv)) = 0;
  Ginv(isinf(Ginv)) = 0;

  Uff = Ginv.*Yd;
  uff = ifft(ifftshift(Uff));
  uff = repmat(uff,2,1);

  set(h_title,'str',['Iteration = ', num2str(i)]);
  set(h_y,'YData',y);
  set(h_err,'XData',tt(si:end),'YData',yd_truc(si:end) - y(si:end));
  set(h_uff,'YData',uff)
  set(h_Ginv, 'YData',abs(Ginv));
  set(h_Uff,'YData',abs(Uff));
  set(h_er_max,'XData',1:i,'YData',Emax_vector(1:i));
  set(h_er_rms,'XData',1:i,'YData',Erms_vector(1:i));
  pause(0.0001)
end

end
