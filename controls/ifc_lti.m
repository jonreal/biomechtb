function [Emax,Erms] = ifc_lti(sys, tt, yd, snr,  ...
                                cutOffFreq, resestIc, numOfIter)
% Iterative feedforward control (ifc)
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
Emax =@(yd, y) norm(yd - y,inf);
Erms =@(yd, y) norm(yd - y,2);

% Frequency vector
f = (0:LL/2)*fs/LL;
f = f(:);

%rho = 0*f + 0.5*(1 + 1j);
rho = 0*f + 0.1;

rho = rho.*(abs(f) < cutOffFreq);

% Start index for fft (half of the yd)
si = LL + 1;

% FFT fo desired traj.
Yd = fft(yd(si:end));
Yd = Yd(1:LL/2 + 1);
Yd(2:end) = Yd(2:end);

%Yd = Yd.*(abs(f) < cutOffFreq);
%  figure;
%    plot(f,abs(Yd))
%    pause

% Find peaks
[~,index] = findpeaks(abs(Yd),'SORTSTR','descend','npeaks',20);

% Truncated
yd_truc = ifft([Yd(1); Yd(2:end); conj(flipud(Yd(2:end-1)))]);
yd_truc = yd_truc(:);
yd_truc = repmat(yd_truc,2,1);

%figure; hold all;
%  plot(tt,yd_truc,'r')
%  plot(tt,yd,'--k')
%  pause

% Bode
frange = logspace(-4,4,1000);
[mag, phi] = bode(sys,frange*2*pi);
mag = mag(:);
phi = phi(:);

fq = abs(f(index));
[mag_yd, phi_yd] = bode(sys,fq*2*pi:);
mag_yd = mag_yd(:);
phi_yd = phi_yd(:);

% ---- Simulation ---- %

Uff = 0.0*Yd;
uff = 0.0.*Uff(1:end-1);
uff = repmat(uff,4,1);

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
E_prev = f*0 + inf;
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
  E = Yd - Yk;

  for j=1:(LL/2 + 1)
    if abs(E(j)) > abs(E_prev(j))
%      phi = angle(E(j)) - angle(E_prev(j));
%      xx = cos(phi);
%      yy = sin(phi);
%      rho(j) = (-yy + 1j*xx)*rho(j)*0.9;
%
%      rho(j) = 1j*rho(j)*0.5;

      rho(j) = rho(j)*0.5;

      Uff(j) = Uff_prev(j);

   % elseif abs(E(j)) < abs(E_prev(j))
    %  rho(j) = rho(j)*1.25;
    end
  end
  E_prev = E;
  Uff_prev = Uff;

  Uff = Uff + rho.*E;

  uff = ifft([Uff(1); Uff(2:end); conj(flipud(Uff(2:end-1)))]);
  uff = repmat(uff,2,1);

  set(h_title,'str',['Iteration = ', num2str(i)]);
  set(h_y,'YData',y);
  set(h_err,'XData',tt(si:end),'YData',yd_truc(si:end) - y(si:end));
  set(h_uff,'YData',uff)
  set(h_E, 'YData',abs(E));
  set(h_Uff,'Ydata',abs(Uff));
  set(h_rho,'YData',abs(rho));
  set(h_er_max,'XData',1:i,'YData',Emax_vector(1:i));
  set(h_er_rms,'XData',1:i,'YData',Erms_vector(1:i));
  pause(0.0001)
end

Emax = Emax_vector(:);
Erms = Erms_vector(:);

end




