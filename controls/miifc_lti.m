function [Emax, Erms] = miifc_lti(sys, tt, yd, snr,  ...
                                cutOffFreq, resestIc, numOfIter)
% Model free inversion based iterative control
% Kim and Zou, 2013
%
% function [uff]  = ifc_lti(sys, tt, yd, snr,  ...
%                                cutOffFreq, resestIc, numOfIter)
%
%
[~,n] = size(sys.a);

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
f = (0:LL/2)*fs/LL;
f = f(:);


% Start index for fft (half of the yd)
si = LL + 1;

% FFT fo desired traj.
Yd = fft(yd(si:end));
Yd = Yd(1:LL/2 + 1);

%Yd = Yd.*(abs(f) < cutOffFreq);
%  figure;
%    plot(f,abs(Yd))
%    pause

% Find peaks
[~,index] = findpeaks(abs(Yd),'SORTSTR','descend','npeaks',20);

% Truncated
yd_truc = ifft([Yd(1:end); conj(flipud(Yd(2:end-1)))]);
yd_truc = yd_truc(:);
yd_truc = repmat(yd_truc,2,1);

figure; hold all;
  plot(tt,yd,'ok')
  plot(tt,yd_truc,'+r')
figure;
  plot(tt,yd-yd_truc);
  pause

% Bode
frange = logspace(-4,4,1000);
[mag, phi] = bode(sys,frange*2*pi);
mag = mag(:);
phi = phi(:);

fq = abs(f(index));
[mag_yd, phi_yd] = bode(sys,fq*2*pi);
mag_yd = mag_yd(:);
phi_yd = phi_yd(:);

% ---- Simulation ---- %

uff = 0.00001*yd_truc;
Uff = 0.00001*Yd;


Emax_vector = zeros(1,numOfIter);
Erms_vector = zeros(1,numOfIter);

figure;

  subplot(4,3,1); hold all; semilogx(frange,mag2db(mag),'k');
    semilogx(fq,mag2db(mag_yd),'or');
    plot([cutOffFreq, cutOffFreq], ylim, '--k')
    set(gca,'xscale','log');
    xlim([min(frange),max(frange)]);
    ylabel('Mag (dB)','fontsize',17,'interpreter','latex');
    grid on
    box on

  subplot(4,3,4); hold all;
    semilogx(frange,phi,'k');
    semilogx(fq,phi_yd,'or');
    plot([cutOffFreq, cutOffFreq], ylim, '--k')
    set(gca,'xscale','log');
    xlim([min(frange),max(frange)]);
    ylabel('Phase (deg)','fontsize',17,'interpreter','latex');
    xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
    grid on
    box on

   subplot(4,3,2:3), hold all;
    h_title = title('Iteration = ','fontsize',20);
    set(h_title,'interpreter','latex');
    plot(tt,yd_truc,'r');
    h_y = plot(tt,zeros(size(tt)));
    h_err = plot(tt,zeros(size(tt)));
    ylim([-2,2]);
    xlabel('Time (s)','fontsize',17,'interpreter','latex');
    legend('Desired','Response','Error');
    grid on
    box on

  subplot(4,3,5:6); hold all;
    h_uff = plot(tt,uff,'m');
    xlabel('Time (s)','fontsize',17,'interpreter','latex');
    grid on
    box on

  subplot(4,3,8); hold all;
    h_E = plot(f(1:LL/2+1),0*f(1:LL/2+1));
    xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
    ylabel('$|E|$','fontsize',17,'interpreter','latex');
    grid on
    box on

  subplot(4,3,9); hold all;
    h_Uff = plot(f(1:LL/2+1),0*f(1:LL/2+1));
    xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
    ylabel('$|U_{ff}|$','fontsize',17,'interpreter','latex')
    grid on
    box on

  subplot(4,3,7); hold all;
    h_rho = plot(f(1:LL/2+1),0*f(1:LL/2+1));
    xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
    ylabel('$\rho$','fontsize',17,'interpreter','latex')
    grid on
    box on

  subplot(4,3,11:12); hold all
     h_er_max = plot(0,0,'k');
     h_er_rms = plot(0,0,'r');
     xlabel('Iteration','fontsize',17,'interpreter','latex');
     ylabel('Error','fontsize',17,'interpreter','latex')
     grid on
     box on

  set(gcf,'Position',[20,20,1000,1000]);


x0 = zeros(n,1);
for i=1:numOfIter

  tic
  [y,~,x] = lsim(sys,uff,tt,x0);
  y = awgn(y,snr,'measured');
  toc

  if resestIc
    x0 = zeros(n,1);
  else
    x0 = x(end,:);
  end

  Erms_vector(i) = Erms(yd_truc(si:end),y(si:end));
  Emax_vector(i) = Emax(yd_truc(si:end),y(si:end));

  % Update
  Yk = fft(y(si:end));
  Yk = Yk(1:LL/2 + 1);

  Ginv = Uff./Yk;

  Uff = Ginv.*Yd;
  Uff(isnan(Uff)) = 0;
  Uff(isinf(Uff)) = 0;


  uff = ifft([Uff; conj(flipud(Uff(2:end-1)))]);
  uff = repmat(uff,2,1);


  set(h_title,'str',['Iteration = ', num2str(i)]);
  set(h_y,'YData',y);
  set(h_err,'XData',tt(si:end),'YData',yd_truc(si:end) - y(si:end));
  set(h_uff,'YData',uff)
  set(h_E, 'YData',abs(Ginv));
  set(h_Uff,'Ydata',abs(Uff));
%  set(h_rho,'YData',abs(rho));
  set(h_er_max,'XData',1:i,'YData',Emax_vector(1:i));
  set(h_er_rms,'XData',1:i,'YData',Erms_vector(1:i));
  pause(0.0001)
end

Emax = Emax_vector(:);
Erms = Erms_vector(:);

end




