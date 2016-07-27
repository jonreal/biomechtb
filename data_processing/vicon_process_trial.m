function rtn = vicon_process_trial(trialName,varargin)

  nVarArgs = length(varargin);

  % Defaults
  stackmethod = 'vicon';
  processViconEvents = 1;
  processEmbedded = 1;

  for i=1:2:nVarArgs
    switch varargin{i}
      case 'stackmethod'
        stackmethod = varargin{i+1};
      case 'processViconEvents'
        processViconEvents = varargin{i+1};
      case 'processEmbedded'
        processEmbedded = varargin{i+1};
      otherwise
        fprintf('\n%s option not found!\n',varargin{i});
        return
    end
  end

  % Process inverse model first
  rtn = vicon_process_model(trialName);

  % Process event data, (must pass time vector) if options are set for that
  if (strcmp(stackmethod,'vicon') || (processViconEvents))
    rtn.segmentedGaitCycles = vicon_process_events(trialName,rtn.time);
  end

  % Process subject params
  rtn.subjectParams = vicon_process_subject(trialName);

  % Process embedded data
  if (strcmp(stackmethod,'embedded') || (processEmbedded))
    rtn.emb = embedded_process_data(trialName);
  end

  % Gait cycle
  gaitCycle = linspace(0,100,1001);
  rtn.stats.gaitCycle = gaitCycle;

  % Store stacking method
  rtn.stats.method = stackmethod;

  % Field names
  q_name = fieldnames(rtn.id);
  foot = {'l','r'};

  % Seperate each gait cycle
  for k=1:numel(foot)

    % Heelstike
    if strcmp(stackmethod,'embedded')
      hs_time = rtn.emb.segmentedGaitCycles.(foot{k}).time;
    elseif strcmp(stackmethod,'vicon')
      hs_time = rtn.segmentedGaitCycles.(foot{k}).time;
    end

    for i=1:numel(hs_time(:,1))

      % If recorded hs occurs before id data, skip
      if hs_time(i,1) < rtn.time(1)
        continue;
      end

      % Find time stamp corresponding to HS
      [~,ii] = min(abs(rtn.time - hs_time(i,1)));
      hs1 = ii;

      [~,ii] = min(abs(rtn.time - hs_time(i,2)));
      hs2 = ii;

      rawTime = rtn.time(hs1:hs2);
      rawGaitCylce = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;

      rtn.stats.normalized.(foot{k}).time{i} = rawTime;

      % Interpolate inverse dynamics
      for j=1:numel(q_name)
        rtn.stats.normalized.(foot{k}).(q_name{j}).X(:,i) = ...
          interp1(rawGaitCylce, rtn.id.(q_name{j}).X(hs1:hs2), ...
                  gaitCycle, 'spline');
        rtn.stats.normalized.(foot{k}).(q_name{j}).Y(:,i) = ...
           interp1(rawGaitCylce, rtn.id.(q_name{j}).Y(hs1:hs2), ...
                   gaitCycle, 'spline');
        rtn.stats.normalized.(foot{k}).(q_name{j}).Z(:,i) = ...
           interp1(rawGaitCylce, rtn.id.(q_name{j}).Z(hs1:hs2), ...
                   gaitCycle, 'spline');
      end
    end

    % Mean, std
    for j=1:numel(q_name)
     jointVarMean_X = mean(rtn.stats.normalized.(foot{k}).(q_name{j}).X, 2);
     jointVarMean_Y = mean(rtn.stats.normalized.(foot{k}).(q_name{j}).Y, 2);
     jointVarMean_Z = mean(rtn.stats.normalized.(foot{k}).(q_name{j}).Z, 2);

     jointVarStd_X = std(rtn.stats.normalized.(foot{k}).(q_name{j}).X')';
     jointVarStd_Y = std(rtn.stats.normalized.(foot{k}).(q_name{j}).Y')';
     jointVarStd_Z = std(rtn.stats.normalized.(foot{k}).(q_name{j}).Z')';

     rtn.stats.(foot{k}).(q_name{j}).X = [jointVarMean_X + jointVarStd_X, ...
                                          jointVarMean_X, ...
                                          jointVarMean_X - jointVarStd_X];
     rtn.stats.(foot{k}).(q_name{j}).Y = [jointVarMean_Y + jointVarStd_Y, ...
                                          jointVarMean_Y, ...
                                          jointVarMean_Y - jointVarStd_Y];
     rtn.stats.(foot{k}).(q_name{j}).Z = [jointVarMean_Z + jointVarStd_Z, ...
                                          jointVarMean_Z, ...
                                          jointVarMean_Z - jointVarStd_Z];
   end
 end

  figure;
    subplot(311); hold all;
      plot_std(rtn.stats.gaitCycle,rtn.stats.r.RAnkleAngles.X,0.1*[1 1 1]);
      plot_std(rtn.stats.gaitCycle,rtn.stats.l.LAnkleAngles.X,0.5*[1 1 1]);
      grid on
    subplot(312); hold all;
      plot_std(rtn.stats.gaitCycle,rtn.stats.r.RAnkleMoment.X,0.1*[1 1 1]);
      plot_std(rtn.stats.gaitCycle,rtn.stats.l.LAnkleMoment.X,0.5*[1 1 1]);
      grid on
    subplot(313); hold all;
      plot_std(rtn.stats.gaitCycle,rtn.stats.r.RAnklePower.X,0.1*[1 1 1]);
      plot_std(rtn.stats.gaitCycle,rtn.stats.l.LAnklePower.X,0.5*[1 1 1]);
      grid on
end
