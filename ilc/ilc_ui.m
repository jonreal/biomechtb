function rtn = ilc_ui(varargin)

  nVarArgs = length(varargin);

  % Defauts Settings
  gain = 1;
  stackmethod = 'vicon';
  maxharmonic = 10;
  torque2current = 0.4;
  earlysmooth = 2;
  latesmooth = 20;
  mass = 85;
  maxcurrent = 20;
  isBioRight = 1;

  for i=1:2:nVarArgs
    switch varargin{i}
      case 'isBioRight'
        isBioRight = varargin{i+1};
      case 'maxcurrent'
        maxcurrent = varargin{i+1}
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
      case 'earlysmooth'
        earlysmooth = varargin{i+1};
      case 'latesmooth'
        latesmooth = varargin{i+1};
      otherwise
        fprintf('\n%s option not found!\n',varargin{i})
        rtn = [];
        return
    end
  end

  % S_0 struct
  S_0.U_k = zeros(501,1);
  S_0.U_k_m_1 = zeros(501,1);
  S_0.E_k_m_1 = ones(501,1).*inf;
  S_0.rho_k_m_1 = ones(501,1);

  % Apply harmonic cutoff -> maxharmonic + 2
  S_0.rho_k_m_1((maxharmonic+2):end) = 0;
  rtn.S{1} = S_0;

  % Filter for uff
  [b,a] = butter(1,(2*10)/1000);

  % Time domain signal length
  L = 1001;

  h_torqe = figure;
  h_error_frq = figure;
  h_rho = figure;
  h_error_trial = figure;
  h_torque_only = figure;


  fprintf('\n----------------------\n')
  fprintf('Iterative Learning UI\n')
  fprintf('-----------------------\n')
  fprintf('\nSettings:\n')
  fprintf('\tBio is on the right side = %i\n',isBioRight);
  fprintf('\tStack Method = %s\n',stackmethod)
  fprintf('\tMaximum harmonic = %i\n',maxharmonic)
  fprintf('\tAnkle torque to motor current = %.3f\n',torque2current)
  fprintf('\tLearning Gain = %.3f\n',gain)
  fprintf('\tEarly smoothing = %i%%\n', earlysmooth)
  fprintf('\tLate smoothing = %i%%\n', latesmooth)

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

    while(1)

      % Iterative Learning
      fprintf('\n\tLearning...\n');
      rtn.S{k+1} = iterativeLearning(rtn.T{k},rtn.S{k},gain,isBioRight);

      % Filter uff
      u_k_filt = rtn.S{k+1}.u_k;
      u_k_filt(1:(round(earlysmooth*10))) = 0;
      u_k_filt(numel(u_k_filt) - (round(latesmooth*10)):end) = 0;
      u_k_filt = filtfilt(b,a,u_k_filt);

      u1 = mass*torque2current*rtn.S{k+1}.u_k;
      u2 = mass*torque2current*u_k_filt;
      if sum(abs(u2) > maxcurrent) > 0
        fprintf('\n---Motor will saturate!---\n');
        u2(u2 > maxcurrent) = maxcurrent;
        u2(u2 < -maxcurrent) = -maxcurrent;
      end

      figure(h_torque_only); clf; hold all;
        plot(rtn.S{k+1}.gaitCycle, u1,'k');
        plot(rtn.S{k+1}.gaitCycle, u2,'r');
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
              earlysmooth = usr_input;
            case 3
              fprintf('\t\t\t');
              usr_input = input('Enter new latesmooth: ');
              latesmooth = usr_input;
            case 4
              fprintf('\tLearning Gain = %.3f\n',gain)
              fprintf('\tEarly smoothing = %i%%\n', earlysmooth)
              fprintf('\tLate smoothing = %i%%\n', latesmooth)
              break;
            otherwise
              warning('Not an option');
          end
        end
      else
        break;
      end
    end

    rtn.S{k+1}.u_k_filt = u_k_filt;
    rtn.S{k+1}.u_k_lut = u2;

    % Torque
    figure(h_torqe); clf; hold all;
      for i=1:k;
        plot3(rtn.S{i+1}.gaitCycle.*0 + (i-1), ...
              rtn.S{i+1}.gaitCycle, ...
              rtn.S{i+1}.yd_k_m_1,'k')
        plot3(rtn.S{i+1}.gaitCycle.*0 + (i-1), ...
              rtn.S{i+1}.gaitCycle, ...
              rtn.S{i+1}.y_k_m_1,'r')
        plot3(rtn.S{i+1}.gaitCycle.*0 + (i-1), ...
              rtn.S{i+1}.gaitCycle, ...
              rtn.S{i+1}.u_k_filt,'g')
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
        stem3(rtn.S{i+1}.f(1:11).*0 + (i-1), ...
              rtn.S{i+1}.f(1:11), ...
              abs(rtn.S{i+1}.Yd_k_m_1(1:11)),'k')
        stem3(rtn.S{i+1}.f(1:11).*0 + (i-1), ...
              rtn.S{i+1}.f(1:11), ...
              abs(rtn.S{i+1}.Y_k_m_1(1:11)),'r')
        stem3(rtn.S{i+1}.f(1:11).*0 + (i-1), ...
              rtn.S{i+1}.f(1:11), ...
              abs(rtn.S{i+1}.E_k_m_1(1:11)),'g')
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
        stem3(rtn.S{i+1}.f(1:11).*0 + (i-1), ...
              rtn.S{i+1}.f(1:11), ...
              abs(rtn.S{i+1}.rho_k_m_1(1:11)),'k')
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
        stem((i-1),rtn.S{i+1}.e_1,'k')
        stem((i-1),rtn.S{i+1}.e_2,'r')
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
        fprintf(fid,'%f\n',rtn.S{k+1}.u_k_lut);
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
