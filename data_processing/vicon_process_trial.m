function rtn = vicon_process_trial(trialName,varargin)

  nVarArgs = length(varargin);

  % Defaults
  stackmethod = 'embedded';
  processViconEvents = 0;
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
    temp = vicon_process_events(trialName,rtn.time);
    rtn.segmentedGaitCycles = temp.segmentedGaitCycles;
    rtn.gaitEvents = temp.gaitEvents; 
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

  % Embedded system feild names (we don't need to interpolate all of them)
  if (strcmp(stackmethod,'embedded') || (processEmbedded))
    q_name_emb = fieldnames(rtn.emb.data);
    q_name_emb = {q_name_emb{12:end-6}}';
  end

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

      % Find time stamp corresponding to HS (vicon data)
      [~,ii] = min(abs(rtn.time - hs_time(i,1)));
      hs1 = ii;

      [~,ii] = min(abs(rtn.time - hs_time(i,2)));
      hs2 = ii;

      rawTime = rtn.time(hs1:hs2);
      rawGaitCylce = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;

      rtn.stats.normalized.(foot{k}).time{i} = rawTime;


      if (strcmp(stackmethod,'embedded') || (processEmbedded))
        % Find time stamp corresponding to HS (embedded data)
        [~,ii] = min(abs(rtn.emb.data.time - hs_time(i,1)));
        hs1_emb = ii;

        [~,ii] = min(abs(rtn.emb.data.time - hs_time(i,2)));
        hs2_emb = ii;

        rawTime = rtn.emb.data.time(hs1_emb:hs2_emb);
        rawGaitCylce_emb = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;

        rtn.emb.stats.normalized.(foot{k}).time{i} = rawTime;
      end

      % Interpolate inverse dynamics (vicon)
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

      % Interpolate embedded system
      if (strcmp(stackmethod,'embedded') || (processEmbedded))
        for j=1:numel(q_name_emb)
          rtn.emb.stats.normalized.(foot{k}).(q_name_emb{j})(:,i) = ...
            interp1(rawGaitCylce_emb, ...
                  rtn.emb.data.(q_name_emb{j})(hs1_emb:hs2_emb), ...
                  gaitCycle, 'spline');
        end
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


    if (strcmp(stackmethod,'embedded') || (processEmbedded))
     for j=1:numel(q_name_emb)
       varMean = mean(rtn.emb.stats.normalized.(foot{k}).(q_name_emb{j}), 2);
       varStd = std(rtn.emb.stats.normalized.(foot{k}).(q_name_emb{j})')';

       rtn.emb.stats.(foot{k}).(q_name_emb{j}) = [varMean + varStd, ...
                                                  varMean, ...
                                                  varMean - varStd];
     end
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
