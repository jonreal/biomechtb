function rtn = ilc_poweredAnkle_tui(varargin)
% ilc_poweredAnkle_tui  Text based user interface for iterative learning
%                       control.
%
%   ilc_poweredAnkle_tui() default configuration.
%
%   ilc_poweredAnkle_tui(Name,Value, ...) Name, Value pairs for adjustable
%   settings.
%
%   Parameters: (defaults)
%
%   gain              -  (1), learning gain, value is positive numerical value
%   stackmethod       -  (embedded), averaging source,
%                         value can be [vicon | embedded]
%   maxharmonic       -  (10), max harmonic to learn, positive interger
%   torque2current    -  (0.4), torque to current confersion factor, numerical value
%                         this value is the inverse plant model
%   smoothing         -  ([95,5])traj. smoothing, value is 2-element vector [b,a]
%                         b - % gait to start smoothing
%                         a - % gait to stop smoothing
%                         b > a, (ie, b=90%, a=10%, smooth from 90% - 10% gait)
%   mass              -  (85) subject mass
%   maxcurrent        -  (20) saturation current, numerical value
%   prostheticSide    -  (left) side with prosthetic
%
  nVarArgs = length(varargin);

  % Defauts Settings
  gain = 1;
  stackmethod = 'embedded';
  maxharmonic = 10;
  torque2current = 0.4;
  smooth_vector = [90,10];
  mass = 85;
  maxcurrent = 20;
  prostheticSide = 'left';
  gaitCycle = 0:0.1:(100 - 0.1);

  for i=1:2:nVarArgs
    switch varargin{i}
      case 'prostheticSide'
        prostheticSide = varargin{i+1};
      case 'maxcurrent'
        maxcurrent = varargin{i+1};
      case 'mass'
        mass = varargin{i+1};
      case 'stackmethod'
        stackmethod = varargin{i+1};
      case 'maxharmonic'
        maxharmonic = varargin{i+1};
      case 'torque2current'
        torque2current = varargin{i+1};
      case 'gain'
        gain = varargin{i+1};
      case 'smoothing'
        smooth_vector = varargin{i+1};
      otherwise
        fprintf('\n%s option not found!\n',varargin{i})
        rtn = [];
        return
    end
  end

  % Filter for uff smoothing
  [b,a] = butter(1,(2*20)/1000);

  %----------------------------------------------------------------------------
  % Plots
  % ---------------------------------------------------------------------------
  figure;

  % |E_bar__k| v freq.
    subplot(3,4,5); hold all;
      h_E_bar_v_f_p_gamma = stem(0:maxharmonic, 0.*(0:maxharmonic),'--r');
      h_E_bar_v_f = stem(0:maxharmonic, 0.*(0:maxharmonic),'k');
      ylabel('$\| \bar{E}_k \|$','fontsize',17,'interpreter','latex');
      xlabel('Harmonic','fontsize',17,'interpreter','latex');
      grid on; box on;

    % |U_bar_k| v freq.
    subplot(3,4,9), hold all;
      h_U_bar_v_f = stem(0:maxharmonic, 0.*(0:maxharmonic),'k');
      ylabel('$\| \bar{U}_k \|$','fontsize',17,'interpreter','latex');
      xlabel('Harmonic','fontsize',17,'interpreter','latex');
      grid on; box on;

    % rho v freq
    subplot(3,4,2); hold all;
      h_rho_v_f = stem(0:maxharmonic, 0.*(0:maxharmonic),'k');
      ylabel('$\rho$','fontsize',17,'interpreter','latex');
      xlabel('Harmonic','fontsize',17,'interpreter','latex');
      grid on; box on;

    % |E_k| v freq.
    subplot(3,4,6); hold all;
      h_E_v_f = stem(0:maxharmonic, 0.*(0:maxharmonic),'k');
      ylabel('$\| E_k \|$','fontsize',17,'interpreter','latex');
      xlabel('Harmonic','fontsize',17,'interpreter','latex');
      grid on; box on;

    % |U_k| v freq.
    subplot(3,4,10), hold all;
      h_U_v_f = stem(0:maxharmonic, 0.*(0:maxharmonic),'k');
      ylabel('$\| U_k \|$','fontsize',17,'interpreter','latex');
      xlabel('Harmonic','fontsize',17,'interpreter','latex');
      grid on; box on;

    % y v t
    subplot(3,4,3:4); hold all;
      h_title = title(['Iteration = ', num2str(0)], ...
                       'interpreter','latex','fontsize',17);
      h_y_v_t = plot(gaitCycle,0.*gaitCycle,'k');
      h_yd_v_t = plot(gaitCycle,0.*gaitCycle,'r');
      h_u_k_v_t = plot(gaitCycle,0.*gaitCycle,'g');
      ylabel('(Nm)','fontsize',17,'interpreter','latex');
      xlabel('\% gait','fontsize',17,'interpreter','latex');
      grid on; box on

    % u v t
    subplot(3,4,7:8); hold all;
      h_u_v_t = plot(gaitCycle,0.*gaitCycle,'m');
      h_u_filt_v_t = plot(gaitCycle,0.*gaitCycle,'g');
      xlabel('\% gait','fontsize',17,'interpreter','latex');
      ylabel('(A)','fontsize',17,'interpreter','latex');
      grid on; box on;

    % e v t | f
    subplot(3,4,11:12); hold all
       h_er_inf_t = plot(0,0,'k');
       h_er_2_t = plot(0,0,'r');
       h_er_inf_f = plot(0,0,'--b');
       h_er_2_f = plot(0,0,'--g');
       h_er_inf_t_min = plot(0,0,'ok');
       h_er_2_t_min = plot(0,0,'or');
       h_er_inf_f_min = plot(0,0,'ob');
       h_er_2_f_min = plot(0,0,'og');
       xlabel('Iteration','fontsize',17,'interpreter','latex');
       ylabel('Error','fontsize',17,'interpreter','latex')
       grid on; box on;
       l = legend('$\|e\|_\infty$','$\|e\|_2$','$\|E\|_\infty$','$\|E\|_2$');
       set(l,'interpreter','latex','box','off');

    set(gcf,'Position',[20,20,1100,1100]);

  fprintf('\n----------------------\n')
  fprintf('Iterative Learning UI\n')
  fprintf('-----------------------\n')
  fprintf('\nSettings:\n')
  fprintf('\tProsthetic Side = %s\n',prostheticSide);
  fprintf('\tStack Method = %s\n',stackmethod)
  fprintf('\tMaximum harmonic = %i\n',maxharmonic)
  fprintf('\tAnkle torque to motor current = %.3f\n',torque2current)
  fprintf('\tMass = %i\n',mass)
  fprintf('\tLearning Gain = %.3f\n',gain)
  fprintf('\tSmoothing = %i%% - %i%%\n', smooth_vector)

  % Error Vectors
  Einf_t = [];
  E2_t = [];
  Einf_f = [];
  E2_f = [];

  k=1;
  while(1)
    fprintf('\n................\n')
    fprintf('Iteration %i\n',k - 1);
    fprintf('................\n')

    % Get k trial
    while(1)
      fprintf('\n\t');
      trial = input('Enter trial name: ','s');
      fprintf(['\n\t',trial]);
      usr_input = input(' - Is this correct? [y/n] ','s');
      if strcmp('y',usr_input)
        break;
      end
    end

    % Data Analysis
    fprintf('\n\tParsing data...\n');
    rtn.T{k} = vicon_process_trial(trial,'stackmethod',stackmethod);

    % Desired output and output
    if strcmp(prostheticSide,'left')

      % Means
      yd_k = rtn.T{k}.stats.r.RAnkleMoment.X(:,2);
      y_k = rtn.T{k}.stats.l.LAnkleMoment.X(:,2);

      % all normalized trajectories
      yd_k_all = rtn.T{k}.stats.normalized.r.RAnkleMoment.X;
      y_k_all = rtn.T{k}.stats.normalized.l.LAnkleMoment.X;

    else

      % Means
      yd_k = rtn.T{k}.stats.l.LAnkleMoment.X(:,2);
      y_k = rtn.T{k}.stats.r.RAnkleMoment.X(:,2);

      % all normalized trajectories
      yd_k_all = rtn.T{k}.stats.l.normalized.LAnkleMoment.X;
      y_k_all = rtn.T{k}.stats.r.normalized.RAnkleMoment.X;
    end

    % Remove last point (0% == 100%)
    yd_k = yd_k(1:(end-1));
    y_k = y_k(1:(end-1));

    while(1)

      % If first iteration, initialize learning
      if (k==1)
        fprintf('\n\tInitializing learning structures...\n');

        % First iteration, find y_0 std in freq. domain
        L = numel(y_k_all(1:end-1,1));
        f = (0:(L/2));
        Y_k_all = fft(y_k_all(1:end-1,:));
        Y_k_all = Y_k_all(1:(L/2)+1,:);
        Y_k_abs_mean = mean(abs(Y_k_all)')';
        Y_k_abs_std = std(abs(Y_k_all)')';

        gamma_ = 3.*Y_k_abs_std;

        figure;
          title('Standard deviation of output','fontsize',20);
          plot(f,gamma_,'ok');
          xlabel('Harmonic','fontsize',20);
          ylabel('STD','fontsize',20)
          xlim([0,maxharmonic])

          pause

        rtn.S{k} = adaptiveILC(yd_k,y_k,gamma_,gain,0,'init', maxharmonic);
      else
        fprintf('\n\tLearning...\n');
        rtn.S{k} = adaptiveILC(yd_k,y_k,gamma_,gain,rtn.S{k-1});
      end

      % Store the errors
      Einf_t(k) = rtn.S{k}.e_k_inf;
      E2_t(k) = rtn.S{k}.e_k_2;
      Einf_f(k) = rtn.S{k}.E_k_inf;
      E2_f(k) = rtn.S{k}.E_k_2;

      [~,iinf_t] = min(Einf_t(1:k)./Einf_t(1));
      [~,i2_t] = min(E2_t(1:k)./E2_t(1));
      [~,iinf_f] = min(Einf_f(1:k)./Einf_f(1));
      [~,i2_f] = min(E2_f(1:k)./E2_f(1));


      % Filter uff  [Nm/kg]
      u_kp1_filt = rtn.S{k}.u_kp1;
      u_kp1_filt(1:(round(smooth_vector(2)*10))) = 0;
      u_kp1_filt(round(smooth_vector(1)*10):end) = 0;
      u_kp1_filt = filtfilt(b,a,u_kp1_filt);

      % Convert to [A]
      u_kp1_A = mass*rtn.S{k}.u_kp1*torque2current;
      u_kp1_A_filt = mass*u_kp1_filt*torque2current;

      if sum(abs(u_kp1_A_filt) > maxcurrent) > 0
        fprintf('\n---Motor will saturate!---\n');
        saturate_pos_indx = find(u_kp1_A_filt > maxcurrent);
        saturate_neg_indx = find(u_kp1_A_filt < -maxcurrent);

        u_kp1_A_filt(saturate_pos_indx) = maxcurrent;
        u_kp1_A_filt(saturate_neg_indx) = -maxcurrent;

      end

      % Convert to [Nm]
      u_kp1_Nm_filt = u_kp1_A_filt ./ torque2current;
      u_kp1_Nm = u_kp1_A ./ torque2current;

      rtn.S{k}.u_kp1_Nm_filt = u_kp1_Nm_filt;
      rtn.S{k}.u_kp1_Nm = u_kp1_Nm;

      rtn.S{k}.u_kp1_A_filt = u_kp1_A_filt;
      rtn.S{k}.u_kp1_A = u_kp1_A;

      % update plots
      set(h_title,'string',['Iteration = ', num2str(k-1)]);
      set(h_u_v_t,'YData',u_kp1_A);
      set(h_u_filt_v_t,'YData',u_kp1_A_filt);
      if k > 1
        set(h_u_k_v_t, 'YData',rtn.S{k-1}.u_kp1_Nm_filt);
        set(h_U_v_f,'YData',abs(rtn.S{k-1}.U_kp1(1:(maxharmonic+1))));
      end
      set(h_U_bar_v_f,'YData',abs(rtn.S{k}.U_bar_k(1:(maxharmonic+1))));
      set(h_E_bar_v_f,'YData',abs(rtn.S{k}.E_bar_k(1:(maxharmonic+1))));
      set(h_E_bar_v_f_p_gamma,'YData', ...
        abs(rtn.S{k}.E_bar_k(1:(maxharmonic+1))) ...
          + rtn.S{k}.gamma_k(1:(maxharmonic+1)))
      set(h_y_v_t,'YData', rtn.S{k}.y_k.*mass);
      set(h_yd_v_t,'YData',rtn.S{k}.yd_k.*mass);

      set(h_er_inf_t_min,'XData',iinf_t,'YData',Einf_t(iinf_t)./Einf_t(1));
      set(h_er_2_t_min,'XData',i2_t,'YData',E2_t(i2_t)./E2_t(1));
      set(h_er_inf_f_min,'XData',iinf_f,'YData',Einf_f(iinf_f)./Einf_f(1));
      set(h_er_2_f_min,'XData',i2_f,'YData',E2_f(i2_f)./E2_f(1));

      set(h_er_inf_t,'XData',1:k,'YData',Einf_t(1:k)./Einf_t(1));
      set(h_er_2_t,'XData',1:k,'YData',E2_t(1:k)./E2_t(1));
      set(h_er_inf_f,'XData',1:k,'YData',Einf_f(1:k)./Einf_f(1));
      set(h_er_2_f,'XData',1:k,'YData',E2_f(1:k)./E2_f(1));

      set(h_rho_v_f,'YData',rtn.S{k}.rho_k(1:(maxharmonic+1)).*gain);
      set(h_E_v_f,'YData',abs(rtn.S{k}.E_k(1:(maxharmonic+1))));

      % Promt for signal sufficient
      fprintf('\n\tSignal learned!\n');
      fprintf('\t\tIs this signal sufficient?')
      usr_input = input('[y|n] ', 's');
      if strcmp(usr_input,'n')
        while(1)
          fprintf('\t\tDo you want to change a parameter?\n')
          fprintf('\t\t 1 - gain\n')
          fprintf('\t\t 2 - earlysmooth\n')
          fprintf('\t\t 3 - latesmooth\n');
          fprintf('\t\t 4 - done\n')
          while(1)
            fprintf('\t\t\t');
            usr_input = input('Enter response: ');
            if isscalar(usr_input)
              break;
            end
          end

          switch usr_input
            case 1
              fprintf('\t\t\t');
              usr_input = input('Enter new gain value: ');
              gain = usr_input;
            case 2
              fprintf('\t\t\t');
              usr_input = input('Enter new earlysmooth: ');
              smooth_vector(2) = usr_input;
            case 3
              fprintf('\t\t\t');
              usr_input = input('Enter new latesmooth: ');
              smooth_vector(1) = usr_input;
            case 4
              fprintf('\tLearning Gain = %.3f\n',gain)
              fprintf('\tSmoothing = %i%% - %i%%\n', smooth_vector)
              break;
            otherwise
              warning('Not an option');
          end
        end
      else
        break;
      end
    end

    % End of Iteration Promt
    fprintf('\n\t');
    fprintf('Write feedforward signal to file?');
    usr_input = input('[y|n]','s');
    fprintf('\n');

    if strcmp(usr_input,'y')
      writeFlag = 1;
      file = ['./','uff_',num2str(k)];
      if exist(file,'file') == 2
        fprintf('\n\t\tFile already exist');
        usr_input = input(' overwrite? [y/n]','s');
        if strcmp('y',usr_input)
          writeFlag = 1;
        else
          writeFlag = 0;
          fprintf('\n\tSignal not written to file.\n');
        end
      end
      if writeFlag
        fid = fopen(file,'w');
        fprintf(fid,'%f\n',rtn.S{k}.u_kp1_A_filt);
        fclose(fid);
      end
    end

    fprintf('\n\t');
    usr_input = input('Continue with next k? [y/n]','s');
    if strcmp('n',usr_input)
      fprintf('\t\t')
      usr_input = input('Redo k? [y/n] (n will exit ui)','s');
      if strcmp('y',usr_input)
        continue;
      end
      break;
    end

    % Increment
    k = k + 1;
  end
end
