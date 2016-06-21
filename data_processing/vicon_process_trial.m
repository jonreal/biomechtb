function rtn = vicon_process_trial(trialName,varargin)

  nVarArgs = length(varargin);

  if nVarArgs > 1
    if strcmp(varargin{1},'stackmethod')
      stackmethod = varargin{2};
    end
  else
    stackmethod = 'vicon';
  end

  % Process inverse model first
  rtn = vicon_process_model(trialName);

  % Process event data, (must pass time vector)
  rtn.gaitEvents = vicon_process_events(trialName,rtn.time);

  % Process subject params
  rtn.subjectParams = vicon_process_subject(trialName);

  % Process embedded data
%  rtn.emb = embedded_process_data(trialName);

  % Gait cycle
  gaitCycle = linspace(0,100,1001);
  rtn.stats.gaitCycle = gaitCycle;

  % Store stacking method
  rtn.stats.method = stackmethod;

  % Field names
  q_name = fieldnames(rtn.id);
  foot = fieldnames(rtn.gaitEvents);

  % Seperate each gait cycle
  for k=1:numel(foot)

    % Heelstike
    if strcmp(stackmethod,'embedded')
      hs_time = rtn.emb.gaitEvents.(foot{k}).hs_time;
    elseif strcmp(stackmethod,'vicon')
      hs_time = rtn.gaitEvents.(foot{k}).hs_time;
    end

    for i=2:numel(hs_time)

      % If recorded hs occurs before id data, skip
      if hs_time(i-1) < rtn.time(1)
        continue;
      end

      % Find time stamp corresponding to HS
      [~,ii] = min(abs(rtn.time - hs_time(i-1)));
      hs1 = ii;

      [~,ii] = min(abs(rtn.time - hs_time(i)));
      hs2 = ii;

      rawTime = rtn.time(hs1:hs2);
      rawGaitCylce = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;

      rtn.stats.normalized.(foot{k}).time{i-1} = rawTime;

      % Interpolate inverse dynamics
      for j=1:numel(q_name)
        rtn.stats.normalized.(foot{k}).(q_name{j}).X(:,i-1) = ...
          interp1(rawGaitCylce, rtn.id.(q_name{j}).X(hs1:hs2), ...
                  gaitCycle, 'spline');
        rtn.stats.normalized.(foot{k}).(q_name{j}).Y(:,i-1) = ...
           interp1(rawGaitCylce, rtn.id.(q_name{j}).Y(hs1:hs2), ...
                   gaitCycle, 'spline');
        rtn.stats.normalized.(foot{k}).(q_name{j}).Z(:,i-1) = ...
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
