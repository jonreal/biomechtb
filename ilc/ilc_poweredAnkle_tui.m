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
  [b,a] = butter(1,(2*10)/1000);

  % Figures
  h_torqe = figure;
  h_error_frq = figure;
  h_rho = figure;
  h_error_trial = figure;
  h_torque_only = figure;

  fprintf('\n----------------------\n')
  fprintf('Iterative Learning UI\n')
  fprintf('-----------------------\n')
  fprintf('\nSettings:\n')
  fprintf('\tProsthetic Side = %s\n',prostheticSide);
  fprintf('\tStack Method = %s\n',stackmethod)
  fprintf('\tMaximum harmonic = %i\n',maxharmonic)
  fprintf('\tAnkle torque to motor current = %.3f\n',torque2current)
  fprintf('\tLearning Gain = %.3f\n',gain)
  fprintf('\tSmoothing = %i%% - %i%%\n', smooth_vector)

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
      yd_k = rtn.T{k}.stats.r.RAnkleMoment.X(:,2);
      y_k = rtn.T{k}.stats.l.LAnkleMoment.X(:,2);
    else
      yd_k = rtn.T{k}.stats.l.LAnkleMoment.X(:,2);
      y_k = rtn.T{k}.stats.r.RAnkleMoment.X(:,2);
    end

    while(1)

      % If first iteration, initialize learning
      if (k==1)
        fprintf('\n\tInitializing learning structures...\n');
        rtn.S{k} = adaptiveILC(yd_k,y_k,gain,0, ...
                               'init',[1,maxharmonic]);
      else
        fprintf('\n\tLearning...\n');
        rtn.S{k} = adaptiveILC(yd_k,y_k,gain,rtn.S{k-1});
      end

      % Filter uff
      u_k_filt = rtn.S{k}.u_kp1;
      u_k_filt(1:(round(smooth_vector(2)*10))) = 0;
      u_k_filt(round(smooth_vector(1)*10):end) = 0;
      u_k_filt = filtfilt(b,a,u_k_filt);

      u1 = mass*torque2current*rtn.S{k}.u_kp1;
      u2 = mass*torque2current*u_k_filt;

      if sum(abs(u2) > maxcurrent) > 0
        fprintf('\n---Motor will saturate!---\n');
        u2(u2 > maxcurrent) = maxcurrent;
        u2(u2 < -maxcurrent) = -maxcurrent;
      end

      figure(h_torque_only); clf; hold all;
        plot(gaitCycle,u1,'k');
        plot(gaitCycle,u2,'r');
        title('Learned Torque Signal','fontsize',20);
        xlabel('% gait','fontsize',20);
        ylabel('Current (A)','fontsize',20);

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
    rtn.S{k}.u_kp1_filt = u_k_filt;
    rtn.S{k}.u_kp1_lt = u2;

    % Torque
    figure(h_torqe); clf; hold all;
      for i=1:k;
        plot3(gaitCycle.*0 + (i-1), ...
              gaitCycle, ...
              rtn.S{i}.yd_k,'k')
        plot3(gaitCycle.*0 + (i-1), ...
              gaitCycle, ...
              rtn.S{i}.y_k,'r')
        plot3(gaitCycle.*0 + (i-1), ...
              gaitCycle, ...
              rtn.S{i}.u_kp1_filt,'g')
      end
      grid on
      if(k>1)
        set(gca,'Xtick',[1:k])
        xlim([0,k-1])
      end
      xlabel('Iteration','fontsize',20);
      ylabel('% gait','fontsize',20);
      zlabel('Nm/kg','fontsize',20);
      title('Torque','fontsize',20);
      view(40,30);

    % Yd and Y plots
    figure(h_error_frq); clf; hold all;
      for i=1:k;
        stem3(rtn.S{i}.f(1:maxharmonic+1).*0 + (i-1), ...
              rtn.S{i}.f(1:maxharmonic+1), ...
              abs(rtn.S{i}.Yd_k(1:maxharmonic+1)),'k')
        stem3(rtn.S{i}.f(1:maxharmonic+1).*0 + (i-1), ...
              rtn.S{i}.f(1:maxharmonic+1), ...
              abs(rtn.S{i}.Y_k(1:maxharmonic+1)),'r')
        stem3(rtn.S{i}.f(1:maxharmonic+1).*0 + (i-1), ...
              rtn.S{i}.f(1:maxharmonic+1), ...
              abs(rtn.S{i}.E_k(1:maxharmonic+1)),'g')
      end
      grid on
      ylim([0,10])
      if(k>1)
        set(gca,'Xtick',[1:k])
        xlim([0,k-1])
      end
      xlabel('Iteration','fontsize',20);
      ylabel('f/f_0','fontsize',20);
      zlabel('Mag','fontsize',20);
      title('Error','fontsize',20);
      view(40,30);

    % Rho
    figure(h_rho); clf; hold all;
      for i=1:k;
        stem3(rtn.S{i}.f(1:maxharmonic+1).*0 + (i-1), ...
              rtn.S{i}.f(1:maxharmonic+1), ...
              abs(rtn.S{i}.rho_k(1:maxharmonic+1)),'k')
      end
      grid on
      ylim([0,10])
      if(k>1)
        set(gca,'Xtick',[1:k])
        xlim([0,k-1])
      end
      xlabel('Iteration','fontsize',20);
      ylabel('f/f_0','fontsize',20);
      zlabel('Mag','fontsize',20);
      title('Learning Weights','fontsize',20);
      view(40,30);

    % Error
    figure(h_error_trial); clf; hold all;
      for i=1:k;
        stem((i-1),rtn.S{i}.e_k_1,'k')
        stem((i-1),rtn.S{i}.e_k_2,'r')
      end
      grid on
      if(k>1)
        set(gca,'Xtick',[1:k])
        xlim([0,k-1])
      end
      xlabel('Iteration','fontsize',20);
      ylabel('Error','fontsize',20);
      title('Error','fontsize',20);

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
        fid = fopen(['./','uff_',num2str(k)],'w');
        fprintf(fid,'%f\n',rtn.S{k}.u_kp1_lt);
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
