function rtn = ilc_ex_2

  global stopFlag;
  stopFlag = 0;

  % Input parameters
  fd = 1;
  numOfPer = 2;
  numOfIter = 100000;
  maxHarmonic = 10;
  gain = 0.1;

  % Time vecotor
  fs = 1000*fd;
  tt = 0:(1/fs):(2*numOfPer*(1/fd)-(1/fs));
  tt = tt(:);

  % System Parameters
  fn = 1;
  zeta = 0.1;
  snr = inf;

  % System
  wn = fn*2*pi;
  s = tf('s');
  sys = (wn^2/(s^2 + 2*zeta*wn*s + wn^2));
  sys_ss = ss(sys);

  ww = 0.001*2*pi;

  A =@(t) [sin(ww.*t) 0.1;
            0.1 cos(ww.*t)];
  B =@(t) [cos(2.*ww.*t)];
  C = [1 1];

  dynm =@(t,x,u) A(t)*x + B(t)*u

  % Desired Trajectory
  Td = 1/fd;
  yd = (mod(tt,Td) < 1/6*Td) * 0 ...
    + ((mod(tt,Td) >= 1/6*Td) & (mod(tt,Td) < 2/6*Td)).*(6/Td*mod(tt,Td) - 1) ...
    + ((mod(tt,Td) >= 2/6*Td) & (mod(tt,Td) < 4/6*Td)).*(-6/Td*mod(tt,Td) + 3) ...
    + ((mod(tt,Td) >= 4/6*Td) & (mod(tt,Td) < 5/6*Td)).*(6/Td*mod(tt,Td) - 5);
  LL = numel(yd);
  L = LL/2;
  f = (0:(L/2)).*fs/L;
  f = f(:);




  % Bode
  frange = logspace(-4,4,1000);
  [mag, phi] = bode(sys,frange*2*pi);
  mag = mag(:);
  phi = phi(:);

  fq = fd*[0:maxHarmonic];
  [mag_yd, phi_yd] = bode(sys,fq*2*pi);
  mag_yd = mag_yd(:);
  phi_yd = phi_yd(:);


  % Plots
  figure;

    subplot(4,3,1); hold all;
      semilogx(frange,mag2db(mag),'k');
      semilogx(fq,mag2db(mag_yd),'or');
      plot([fd*maxHarmonic, fd*maxHarmonic], ylim, '--k')
      set(gca,'xscale','log');
      xlim([min(frange),max(frange)]);
      ylabel('Mag (dB)','fontsize',17,'interpreter','latex');
      grid on
      box on

    subplot(4,3,4); hold all;
      semilogx(frange,phi,'k');
      semilogx(fq,phi_yd,'or');
      plot([fd*maxHarmonic, fd*maxHarmonic], ylim, '--k')
      set(gca,'xscale','log');
      xlim([min(frange),max(frange)]);
      ylabel('Phase (deg)','fontsize',17,'interpreter','latex');
      xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
      grid on
      box on

     subplot(4,3,2:3), hold all;
      h_title = title('Iteration = ','fontsize',20);
      set(h_title,'interpreter','latex');
      plot(tt,yd,'r');
      h_y = plot(tt,zeros(size(tt)));
      h_err = plot(tt,zeros(size(tt)));
      ylim([-2,2]);
      xlabel('Time (s)','fontsize',17,'interpreter','latex');
      legend('Desired','Response','Error');
      grid on
      box on

    subplot(4,3,5:6); hold all;
      h_uff = plot(tt,yd.*0,'m');
      xlabel('Time (s)','fontsize',17,'interpreter','latex');
      grid on
      box on

    subplot(4,3,8); hold all;
      h_E = stem(fq,0*fq);
      xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
      ylabel('$|E|$','fontsize',17,'interpreter','latex');
      grid on
      box on

    subplot(4,3,9); hold all;
      h_Uff = stem(fq,0*fq);
      xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
      ylabel('$|U_{ff}|$','fontsize',17,'interpreter','latex')
      grid on
      box on

    subplot(4,3,7); hold all;
      h_rho = stem(fq,0*fq);
      xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
      ylabel('$\rho$','fontsize',17,'interpreter','latex')
      grid on
      box on

    subplot(4,3,11:12); hold all
       h_er_inf_t = plot(0,0,'k');
       h_er_2_t = plot(0,0,'r');
       h_er_inf_f = plot(0,0,'--b');
       h_er_2_f = plot(0,0,'--g');
       h_er_inf_t_min = plot(0,0,'ok');
       h_er_2_t_min = plot(0,0,'or');
       h_er_inf_f_min = plot(0,0,'ok');
       h_er_2_f_min = plot(0,0,'or');
       xlabel('Iteration','fontsize',17,'interpreter','latex');
       ylabel('Error','fontsize',17,'interpreter','latex')
       grid on
       box on
       legend('INF_t','2_t','INF_f','2_f')

    set(gcf,'Position',[20,20,1000,1000]);

    hbutton=uicontrol(gcf,'style','pushbutton',...
                      'string','End',...
                       'callback',@callBack);
  % Simulation
  for i=1:numOfIter

    if (i==1)

      % No need to simulate

      % Learn last half of trajectory
      yd_ = yd(LL/2+1:end);
      yd_stack = reshape(yd_,Td*fs,numOfPer);
      yd_ = mean(yd_stack,2);
      y_ = 0.*yd_;
      yy = yd.*0;
      u = yd.*0;

      %S{i} = adaptiveILC(yd_, y_, gain, 0, 'init', maxHarmonic);
      S{i} = adaptivePhaseILC(yd_, y_, gain, 0, 'init', maxHarmonic);

    else
      % Simulate with kth control signal
      u = S{i-1}.u_kp1;
      u = repmat(u, 2*numOfPer, 1);

      [yy,~,xx] = lsim(sys_ss, u, tt, [0;0]);
      yy = awgn(yy,snr,'measured');

      uu =@(t) interp1(tt,u,t);
      [~,x] = ode45(@(t,x) dynm(t,x,uu(t)), tt, [0;0])
      yy = x(:,1);

      % Learn last half of trajectory
      y_ = yy(LL/2+1:end);
      y_stack = reshape(y_,Td*fs,numOfPer);
      y_ = mean(y_stack,2);

      %S{i} = adaptiveILC(yd_, y_, gain, S{i-1});
      S{i} = adaptivePhaseILC(yd_, y_, gain, S{i-1});

    end

    Einf_t(i) = S{i}.e_k_inf;
    E2_t(i) = S{i}.e_k_2;
    Einf_f(i) = S{i}.E_k_inf;
    E2_f(i) = S{i}.E_k_2;

    [~,iinf_t] = min(Einf_t(1:i)./Einf_t(1));
    [~,i2_t] = min(E2_t(1:i)./E2_t(1));
    [~,iinf_f] = min(Einf_f(1:i)./Einf_f(1));
    [~,i2_f] = min(E2_f(1:i)./E2_f(1));

    % Update Plots
    set(h_title,'str',['Iteration = ', num2str(i)]);
    set(h_y,'YData',yy);
    set(h_err,'XData',tt(L:end),'YData',yd(L:end) - yy(L:end));
    set(h_uff,'YData',u)
    set(h_E, 'YData',abs(S{i}.E_bar_k(1:(maxHarmonic+1))));
    set(h_Uff,'Ydata',abs(S{i}.U_bar_k(1:(maxHarmonic+1))));
    set(h_rho,'YData',abs(S{i}.rho_k(1:(maxHarmonic+1))));

    set(h_er_inf_t,'XData',1:i,'YData',Einf_t(1:i)./Einf_t(1));
    set(h_er_2_t,'XData',1:i,'YData',E2_t(1:i)./E2_t(1));
    set(h_er_inf_f,'XData',1:i,'YData',Einf_f(1:i)./Einf_f(1));
    set(h_er_2_f,'XData',1:i,'YData',E2_f(1:i)./E2_f(1));

    set(h_er_inf_t_min,'XData',iinf_t,'YData',Einf_t(iinf_t)./Einf_t(1));
    set(h_er_2_t_min,'XData',i2_t,'YData',E2_t(i2_t)./E2_t(1));
    set(h_er_inf_f_min,'XData',iinf_f,'YData',Einf_f(iinf_f)./Einf_f(1));
    set(h_er_2_f_min,'XData',i2_f,'YData',E2_f(i2_f)./E2_f(1));


    fprintf('.....................\n');
    Einf_t_ = [Einf_t(i) Einf_t(iinf_t)] ./ Einf_t(1)
    E2_t_ = [E2_t(i) E2_t(i2_t)] ./ E2_t(1)
    Einf_f_ = [Einf_f(i) Einf_f(iinf_f)] ./ E2_f(1)
    E2_f_ = [E2_f(i) E2_f(i2_f)] ./ Einf_f(1)

    pause(0.0001)

    if stopFlag
      rtn = S;
      break;
    end
  end
  rtn = S;


  figure;
    subplot(211)
    plot(0:(i-1),mag2db(E2_f./E2_f(1)),'k');
    ylabel('$\frac{\|E_k\|_2}{\|E_0\|_2}$','fontsize',20,'interpreter','latex');
    xlabel('$k$','interpreter','latex','fontsize',20);
    grid on; box on;
    set(gca,'XScale','log')
end

function callBack(hObject,eventdata,handles)
  global stopFlag;
  stopFlag = 1;
end


