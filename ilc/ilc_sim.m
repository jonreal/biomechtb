function rtn = ilc_ex_2

  global stopFlag;
  global startFlag;

  rng(2);
  % Input parameters
  fd = 1;
  snr = 10;
  numOfPer = 2;
  numOfIter = 2000;
  maxHarmonic = 300;
  gain = 1;
  adapt = 'special';
  model = 'no';
  backstep = 'no';
  switchYd = 0;
  ydi = 1;

  stopFlag = 0;
  startFlag = 0;
  doMovie = 0;
  nFrames = 200;

  if doMovie
    mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
    numOfIter = nFrames;
  end

  % Time vecotor
  fs = 1000*fd;
  tt = 0:(1/fs):(2*numOfPer*(1/fd)-(1/fs));
  tt = tt(:);

  % System Parameters
  fn = 1;
  zeta = 0.1;

  % System
  wn = fn*2*pi;
  s = tf('s');
  sys = (wn^2/(s^2 + 2*zeta*wn*s + wn^2));
  sys_ss = ss(sys);
  inv_ss = ((s^2 + 2*zeta*wn*s + wn^2)/wn^2)

  % Complex System
  omegaz = 10*2*pi;
  omegap1 = 1*2*pi;
  omegap2 = 20*2*pi;
  zetaz = 0.7; 
  zetap1 = 0.2; 
  zetap2 = 0.2;
  z =[1 -2*zetaz*omegaz  omegaz^2];
  p1 =[1 2*zetap1*omegap1  omegap1^2];
  p2 =[1 2*zetap1*omegap2  omegap2^2];
  sys = tf(z,conv(p1,p2));
  k = omegap1^2*omegap2^2/omegaz^2;
  num = z*k; 
  den = conv(p1,p2); 
  sys = tf(num,den);

  % Time stuff
  LL = numel(tt);
  L = LL/2;
  f = (0:(L/2)).*fs/L;
  f = f(:);

  % Desired Trajectory 1
  Td = 1/fd;
  yd1 = (mod(tt,Td) < 1/8*Td) * 0 ...
    + ((mod(tt,Td) >= 1/8*Td) & (mod(tt,Td) < 2/8*Td)).*(8/Td*mod(tt,Td) - 1) ...
    + ((mod(tt,Td) >= 2/8*Td) & (mod(tt,Td) < 4/8*Td)).*(-8/Td*mod(tt,Td) + 3) ...
    + ((mod(tt,Td) >= 4/8*Td) & (mod(tt,Td) < 5/8*Td)).*(8/Td*mod(tt,Td) - 5);

  % Desired Trajectory 2
  Td2 = Td/2;
  Ad = 2*pi/(Td2^2);
  wd = 2*pi*(1/Td2);
  yd2 = (Ad/wd).*tt(1:LL/(numOfPer*4)) ...
       - (Ad/wd/wd)*sin(wd*tt(1:LL/(numOfPer*4)));
  yd2 = repmat([yd2; flipud(yd2)],2*numOfPer,1);

  yd_both = [yd1,yd2];

  % Bode
  frange = logspace(-4,4,1000);
  [mag, phi] = bode(sys,frange*2*pi);
  mag = mag(:);
  phi = phi(:);

  fq = fd*[0:maxHarmonic];
  [mag_yd, phi_yd] = bode(sys,fq*2*pi);
  mag_yd = mag_yd(:);
  phi_yd = phi_yd(:);

  [inv_mag,inv_phi] = bode(inv_ss,frange.*2*pi);
  inv_mag = inv_mag(:);
  inv_phi = inv_phi(:);

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
      h_yd = plot(tt,tt.*0,'r');
      h_y = plot(tt,zeros(size(tt)));
      h_err = plot(tt,zeros(size(tt)));
      ylim([-2,2]);
      xlabel('Time (s)','fontsize',17,'interpreter','latex');
      legend('Desired','Response','Error');
      grid on
      box on

    subplot(4,3,5:6); hold all;
      h_uff = plot(tt,tt.*0,'m');
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

    subplot(4,3,10); hold all;
      h_rho_phase = stem(fq,0*fq);
      xlim([0,max(fq)]);
      xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
      ylabel('$\angle \rho$','fontsize',17,'interpreter','latex')
      grid on
      box on

    subplot(4,3,7); hold all;
      h_rho = stem(fq,0*fq);
      xlim([0,max(fq)]);
      xlabel('Frequency (Hz)','fontsize',17,'interpreter','latex');
      ylabel('$|\rho|$','fontsize',17,'interpreter','latex')
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
    hbutton2=uicontrol(gcf,'style','pushbutton',...
                      'string','Start',...
                      'Position', [400 20 120 20], ...
                       'callback',@callBack_2);



    while (1)
      pause(0.01);
      if startFlag == 1;
        break;
      end
    end

  
  yd = yd_both(:,ydi);    

  % Simulation
  for i=1:numOfIter


    if (i==1)

      % No need to simulate

      % Learn last half of trajectory
  %    ydd = awgn(yd,snr,'measured');

      yd_ = yd(LL/2+1:end);

      yd_stack = reshape(yd_,Td*fs,numOfPer);
      yd_ = mean(yd_stack,2);
      y_ = 0.*yd_;
      yy = yd.*0;
      u = yd.*0;

      gamma_ = y_.*0;
      gamma_ = gamma_(1:(numel(gamma_)/2+1));

      L = numel(yd_);
      Yd = fft(yd_stack);
      Yd = Yd(1:(L/2+1),:);

      epsilon = 3.*std(Yd')';

      S{i} = adaptivePhaseILC(yd_,y_,0, ...
                             'init','gain',gain, ...
                             'epsilon', 0,...
                             'model',model,...
                             'adapt',adapt,...
                             'backstep',backstep,...
                             'maxHarm',maxHarmonic, ...
                             'rho_max',inf, ...
                             'zeta',1.45, ...
                             'alpha',0.5);

      ydd = yd;
    else

      % Switch trajectory
      if (mod(i,25) == 0) && (switchYd)
        yd_both = fliplr(yd_both);
        yd = yd_both(:,1);

        % Reinitialize
        S{i-1}.u_kp1 = S{i-1}.u_kp1.*0;
        S{i-1}.U_kp1 = S{i-1}.U_kp1.*0 +0.00001;
       % S{i-1}.E_bar_k = S{i-1}.E_bar_k.*0 + 100000000;
       %
        S{i-1}.E_bar_k = S{i-1}.U_bar_k.*0 + 0.000001;
      end

      % Simulate with kth control signal
      u = S{i-1}.u_kp1;
      u = repmat(u, 2*numOfPer, 1);

      U = fft(u);
      U = U(1:(LL/2+1));

      [yy,~,xx] = lsim(sys_ss, u, tt, [0;0]);
      yy = awgn(yy,snr,'measured');

      %uu =@(t) interp1(tt,u,t);
      %[~,x] = ode45(@(t,x) dynm(t,x,uu(t)), tt, [0;0])
      %yy = x(:,1);

      % Learn last half of trajectory
      y_ = yy(LL/2+1:end);
      y_stack = reshape(y_,Td*fs,numOfPer);
      y_ = mean(y_stack,2);

    %  ydd = awgn(yd,snr,'measured');
      yd_ = yd(LL/2+1:end);
      yd_stack = reshape(yd_,Td*fs,numOfPer);
      yd_ = mean(yd_stack,2);

      S{i} = adaptivePhaseILC(yd_, y_,S{i-1});
  
      ydd = yd;
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
    set(h_yd,'YData',ydd);
    set(h_err,'XData',tt(LL/2+1:end),'YData',yd(LL/2+1:end) - yy(LL/2+1:end));
    set(h_uff,'YData',u)
    set(h_E, 'YData',abs(S{i}.E_bar_k(1:(maxHarmonic+1))));
    set(h_Uff,'Ydata',abs(S{i}.U_bar_k(1:(maxHarmonic+1))));

    set(h_rho,'YData',abs(S{i}.rho_k(1:(maxHarmonic+1))));
    set(h_rho_phase,'YData',angle(S{i}.rho_k(1:(maxHarmonic+1))));

    set(h_er_inf_t,'XData',1:i,'YData',Einf_t(1:i)./Einf_t(1));
    set(h_er_2_t,'XData',1:i,'YData',E2_t(1:i)./E2_t(1));
    set(h_er_inf_f,'XData',1:i,'YData',Einf_f(1:i)./Einf_f(1));
    set(h_er_2_f,'XData',1:i,'YData',E2_f(1:i)./E2_f(1));

    set(h_er_inf_t_min,'XData',iinf_t,'YData',Einf_t(iinf_t)./Einf_t(1));
    set(h_er_2_t_min,'XData',i2_t,'YData',E2_t(i2_t)./E2_t(1));
    set(h_er_inf_f_min,'XData',iinf_f,'YData',Einf_f(iinf_f)./Einf_f(1));
    set(h_er_2_f_min,'XData',i2_f,'YData',E2_f(i2_f)./E2_f(1));


%    fprintf('.....................\n');
%    Einf_t_ = [Einf_t(i) Einf_t(iinf_t)] ./ Einf_t(1);
%    E2_t_ = [E2_t(i) E2_t(i2_t)] ./ E2_t(1);
%    Einf_f_ = [Einf_f(i) Einf_f(iinf_f)] ./ E2_f(1);
%    E2_f_ = [E2_f(i) E2_f(i2_f)] ./ Einf_f(1);
%
    if doMovie
      mov(i) = getframe(gcf);
    end

    pause(0.00001)
    if stopFlag
      rtn = S;
      break;
    end
  end
  rtn = S;


  figure; hold all;
    title(['model = ', model,', adapt = ', adapt, ', backstep = ',backstep], ...
          'fontsize',20, 'interpreter','latex');
    plot(0:(i-1),mag2db(E2_t./E2_t(1)),'k');
    ylabel('$\frac{\|E_k\|_2}{\|E_0\|_2}$','fontsize',20,'interpreter','latex');
    xlabel('$k$','interpreter','latex','fontsize',20);
    grid on; box on;
  %  set(gca,'XScale','log')

    if doMovie
      movie2avi(mov, 'ilc_sim.avi', ...
                     'compression', 'None', ...
                      'quality', 25, ...
                      'fps',3);
    end
end

function callBack(hObject,eventdata,handles)
  global stopFlag;
  stopFlag = 1;
end

function callBack_2(hObject,eventdata,handles)
  global startFlag;
  startFlag = 1;
end


