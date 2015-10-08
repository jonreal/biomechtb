function S = vicon_process_trial(trialName)

  S = vicon_process_model(trialName);
  S.GRF = vicon_process_grf(trialName);
  S.Events = vicon_process_events(S);

  % Interpolate Pedar
  PRSR = pedar_process_trial(trialName,'W');
  S.Pedar.Raw.time = PRSR.time;
  S.Pedar.Raw.l_Pa = PRSR.l_Pa;
  S.Pedar.Raw.r_Pa = PRSR.r_Pa;

  S.Pedar.Interp.l_Pa = interp1(PRSR.time,PRSR.l_Pa,S.time,'nearest');
  S.Pedar.Interp.r_Pa = interp1(PRSR.time,PRSR.r_Pa,S.time,'nearest');

  S.Pedar.sensorTemplate = PRSR.sensorTemplate;

  gaitCycle = linspace(0,100,1001);
  S.Stats.gaitCycle = gaitCycle;

  % Left
  indx = S.Events.L.HS;
  MDF = S.Events.L.MDF;
  TO = S.Events.L.TO;
  labels = fieldnames(S.Model);
  for i=2:numel(indx)

    strtIndx = indx(i-1);
    endIndx = indx(i);

    currMDF = MDF;
    currMDF(currMDF < strtIndx) = [];
    currMDF(currMDF > endIndx) = [];

    currTO = TO;
    currTO(currTO < strtIndx) = [];
    currTO(currTO > endIndx) = [];

    rawTime = S.time(strtIndx:endIndx);
    rawGaitCylce = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;
    S.Stats.LEvents.time{i-1} = rawTime;

    % --- % Gait for TO and MDF
    S.Stats.LEvents.MDF(i-1) = ...
        (S.time(currMDF) - rawTime(1))/(rawTime(end) - rawTime(1))*100;
    S.Stats.LEvents.TO(i-1) = ...
        (S.time(currTO) - rawTime(1))/(rawTime(end) - rawTime(1))*100;

    % --- Kinematics
    for j=1:numel(labels)
      S.Stats.LEvents.(labels{j}).X(:,i-1) = ...
        interp1(rawGaitCylce, S.Model.(labels{j}).X(strtIndx:endIndx), ...
                gaitCycle, 'spline');
      S.Stats.LEvents.(labels{j}).Y(:,i-1) = ...
        interp1(rawGaitCylce, S.Model.(labels{j}).Y(strtIndx:endIndx), ...
                gaitCycle, 'spline');
      S.Stats.LEvents.(labels{j}).Z(:,i-1) = ...
        interp1(rawGaitCylce, S.Model.(labels{j}).Y(strtIndx:endIndx), ...
                gaitCycle, 'spline');
    end

    % --- Pressure
    S.Stats.LEvents.l_Pa(:,:,i-1) = interp1(rawGaitCylce, ...
                                      S.Pedar.Interp.l_Pa(strtIndx:endIndx,:),...
                                      gaitCycle, 'spline');
    S.Stats.LEvents.r_Pa(:,:,i-1) = interp1(rawGaitCylce, ...
                                      S.Pedar.Interp.r_Pa(strtIndx:endIndx,:), ...
                                      gaitCycle, 'spline');

    % --- GRF
    S.Stats.LEvents.GRF.L.Fx(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.L.Fx(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.LEvents.GRF.L.Fy(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.L.Fy(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.LEvents.GRF.L.Fz(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.L.Fz(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.LEvents.GRF.R.Fx(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.R.Fx(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.LEvents.GRF.R.Fy(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.R.Fy(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.LEvents.GRF.R.Fz(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.R.Fz(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
  end


  % Left mean, std
  for j=1:numel(labels)
    jointVarMean_X = mean(S.Stats.LEvents.(labels{j}).X, 2);
    jointVarMean_Y = mean(S.Stats.LEvents.(labels{j}).Y, 2);
    jointVarMean_Z = mean(S.Stats.LEvents.(labels{j}).Z, 2);

    jointVarStd_X = std(S.Stats.LEvents.(labels{j}).X')';
    jointVarStd_Y = std(S.Stats.LEvents.(labels{j}).Y')';
    jointVarStd_Z = std(S.Stats.LEvents.(labels{j}).Z')';

    S.Stats.LMean_std.(labels{j}).X = [jointVarMean_X + jointVarStd_X, ...
                                       jointVarMean_X, ...
                                       jointVarMean_X - jointVarStd_X];
    S.Stats.LMean_std.(labels{j}).Y = [jointVarMean_Y + jointVarStd_Y, ...
                                       jointVarMean_Y, ...
                                       jointVarMean_Y - jointVarStd_Y];
    S.Stats.LMean_std.(labels{j}).Z = [jointVarMean_Z + jointVarStd_Z, ...
                                       jointVarMean_Z, ...
                                       jointVarMean_Z - jointVarStd_Z];
  end

  % --- Pressure
  mean_l_Pa = mean(S.Stats.LEvents.l_Pa,3);
  std_l_Pa = std(S.Stats.LEvents.l_Pa,0,3);
  S.Stats.LMean_std.l_Pa = [{mean_l_Pa + std_l_Pa}, ...
                           {mean_l_Pa}, ...
                           {mean_l_Pa - std_l_Pa}];
  mean_r_Pa = mean(S.Stats.LEvents.r_Pa,3);
  std_r_Pa = std(S.Stats.LEvents.r_Pa,0,3);
  S.Stats.LMean_std.r_Pa = [{mean_r_Pa + std_r_Pa}, ...
                           {mean_r_Pa}, ...
                           {mean_r_Pa - std_r_Pa}];
  % --- GRF
  meanFx = mean(S.Stats.LEvents.GRF.L.Fx,2);
  stdFx = std(S.Stats.LEvents.GRF.L.Fx,0,2);
  S.Stats.LMean_std.GRF.L.Fx = [meanFx + stdFx, ...
                                meanFx, ...
                                meanFx - stdFx];
  meanFy = mean(S.Stats.LEvents.GRF.L.Fy,2);
  stdFy = std(S.Stats.LEvents.GRF.L.Fy,0,2);
  S.Stats.LMean_std.GRF.L.Fy = [meanFy + stdFy, ...
                                meanFy, ...
                                meanFy - stdFy];
  meanFz = mean(S.Stats.LEvents.GRF.L.Fz,2);
  stdFz = std(S.Stats.LEvents.GRF.L.Fz,0,2);
  S.Stats.LMean_std.GRF.L.Fz = [meanFz + stdFz, ...
                                meanFz, ...
                                meanFz - stdFz];

  meanFx = mean(S.Stats.LEvents.GRF.R.Fx,2);
  stdFx = std(S.Stats.LEvents.GRF.R.Fx,0,2);
  S.Stats.LMean_std.GRF.R.Fx = [meanFx + stdFx, ...
                                meanFx, ...
                                meanFx - stdFx];
  meanFy = mean(S.Stats.LEvents.GRF.R.Fy,2);
  stdFy = std(S.Stats.LEvents.GRF.R.Fy,0,2);
  S.Stats.LMean_std.GRF.R.Fy = [meanFy + stdFy, ...
                                meanFy, ...
                                meanFy - stdFy];
  meanFz = mean(S.Stats.LEvents.GRF.R.Fz,2);
  stdFz = std(S.Stats.LEvents.GRF.R.Fz,0,2);
  S.Stats.LMean_std.GRF.R.Fz = [meanFz + stdFz, ...
                                meanFz, ...
                                meanFz - stdFz];

  meanMDF = mean(S.Stats.LEvents.MDF);
  stdMDF = std(S.Stats.LEvents.MDF);
  S.Stats.LMean_std.MDF = [meanMDF + stdMDF, meanMDF, meanMDF - stdMDF];

  meanTO = mean(S.Stats.LEvents.TO);
  stdTO = std(S.Stats.LEvents.TO);
  S.Stats.LMean_std.TO = [meanTO + stdTO, meanTO, meanTO - stdTO];


  % Right
  indx = S.Events.R.HS;
  MDF = S.Events.R.MDF;
  TO = S.Events.R.TO;
  labels = fieldnames(S.Model);
  for i=2:numel(indx)
    strtIndx = indx(i-1);
    endIndx = indx(i);
    rawTime = S.time(strtIndx:endIndx);
    rawGaitCylce = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;
    S.Stats.REvents.time{i-1} = rawTime;

    currMDF = MDF;
    currMDF(currMDF < strtIndx) = [];
    currMDF(currMDF > endIndx) = [];

    currTO = TO;
    currTO(currTO < strtIndx) = [];
    currTO(currTO > endIndx) = [];

    % --- % Gait for TO and MDF
    S.Stats.REvents.MDF(i-1) = ...
        (S.time(currMDF) - rawTime(1))/(rawTime(end) - rawTime(1))*100;
    S.Stats.REvents.TO(i-1) = ...
        (S.time(currTO) - rawTime(1))/(rawTime(end) - rawTime(1))*100;

    % --- Kinematics
    for j=1:numel(labels)
      S.Stats.REvents.(labels{j}).X(:,i-1) = ...
        interp1(rawGaitCylce, S.Model.(labels{j}).X(strtIndx:endIndx), ...
                gaitCycle, 'spline');
      S.Stats.REvents.(labels{j}).Y(:,i-1) = ...
        interp1(rawGaitCylce, S.Model.(labels{j}).Y(strtIndx:endIndx), ...
                gaitCycle, 'spline');
      S.Stats.REvents.(labels{j}).Z(:,i-1) = ...
        interp1(rawGaitCylce, S.Model.(labels{j}).Y(strtIndx:endIndx), ...
                gaitCycle, 'spline');
    end

    % --- Pressure
    S.Stats.REvents.l_Pa(:,:,i-1) = interp1(rawGaitCylce, ...
                                      S.Pedar.Interp.l_Pa(strtIndx:endIndx,:),...
                                      gaitCycle, 'nearest');
    S.Stats.REvents.r_Pa(:,:,i-1) = interp1(rawGaitCylce, ...
                                      S.Pedar.Interp.r_Pa(strtIndx:endIndx,:), ...
                                      gaitCycle, 'nearest');
    % --- GRF
    S.Stats.REvents.GRF.L.Fx(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.L.Fx(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.REvents.GRF.L.Fy(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.L.Fy(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.REvents.GRF.L.Fz(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.L.Fz(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.REvents.GRF.R.Fx(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.R.Fx(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.REvents.GRF.R.Fy(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.R.Fy(strtIndx:endIndx), ...
                                    gaitCycle,'spline');
    S.Stats.REvents.GRF.R.Fz(:,i-1) = interp1(rawGaitCylce, ...
                                    S.GRF.Decimate.R.Fz(strtIndx:endIndx), ...
                                    gaitCycle,'spline');

  end

  % Right mean, std
  for j=1:numel(labels)
    jointVarMean_X = mean(S.Stats.REvents.(labels{j}).X, 2);
    jointVarMean_Y = mean(S.Stats.REvents.(labels{j}).Y, 2);
    jointVarMean_Z = mean(S.Stats.REvents.(labels{j}).Z, 2);

    jointVarStd_X = std(S.Stats.REvents.(labels{j}).X')';
    jointVarStd_Y = std(S.Stats.REvents.(labels{j}).Y')';
    jointVarStd_Z = std(S.Stats.REvents.(labels{j}).Z')';

    S.Stats.RMean_std.(labels{j}).X = [jointVarMean_X + jointVarStd_X, ...
                                       jointVarMean_X, ...
                                       jointVarMean_X - jointVarStd_X];
    S.Stats.RMean_std.(labels{j}).Y = [jointVarMean_Y + jointVarStd_Y, ...
                                       jointVarMean_Y, ...
                                       jointVarMean_Y - jointVarStd_Y];
    S.Stats.RMean_std.(labels{j}).Z = [jointVarMean_Z + jointVarStd_Z, ...
                                       jointVarMean_Z, ...
                                       jointVarMean_Z - jointVarStd_Z];
  end

  % --- Pressure
  mean_l_Pa = mean(S.Stats.REvents.l_Pa,3);
  std_l_Pa = std(S.Stats.REvents.l_Pa,0,3);
  S.Stats.RMean_std.l_Pa = [{mean_l_Pa + std_l_Pa}, ...
                           {mean_l_Pa}, ...
                           {mean_l_Pa - std_l_Pa}];
  mean_r_Pa = mean(S.Stats.REvents.r_Pa,3);
  std_r_Pa = std(S.Stats.REvents.r_Pa,0,3);
  S.Stats.RMean_std.r_Pa = [{mean_r_Pa + std_r_Pa}, ...
                           {mean_r_Pa}, ...
                           {mean_r_Pa - std_r_Pa}];

  % --- GRF
  meanFx = mean(S.Stats.REvents.GRF.L.Fx,2);
  stdFx = std(S.Stats.REvents.GRF.L.Fx,0,2);
  S.Stats.RMean_std.GRF.L.Fx = [meanFx + stdFx, ...
                                meanFx, ...
                                meanFx - stdFx];
  meanFy = mean(S.Stats.REvents.GRF.L.Fy,2);
  stdFy = std(S.Stats.REvents.GRF.L.Fy,0,2);
  S.Stats.RMean_std.GRF.L.Fy = [meanFy + stdFy, ...
                                meanFy, ...
                                meanFy - stdFy];
  meanFz = mean(S.Stats.REvents.GRF.L.Fz,2);
  stdFz = std(S.Stats.REvents.GRF.L.Fz,0,2);
  S.Stats.RMean_std.GRF.L.Fz = [meanFz + stdFz, ...
                                meanFz, ...
                                meanFz - stdFz];

  meanFx = mean(S.Stats.REvents.GRF.R.Fx,2);
  stdFx = std(S.Stats.REvents.GRF.R.Fx,0,2);
  S.Stats.RMean_std.GRF.R.Fx = [meanFx + stdFx, ...
                                meanFx, ...
                                meanFx - stdFx];
  meanFy = mean(S.Stats.REvents.GRF.R.Fy,2);
  stdFy = std(S.Stats.REvents.GRF.R.Fy,0,2);
  S.Stats.RMean_std.GRF.R.Fy = [meanFy + stdFy, ...
                                meanFy, ...
                                meanFy - stdFy];
  meanFz = mean(S.Stats.REvents.GRF.R.Fz,2);
  stdFz = std(S.Stats.REvents.GRF.R.Fz,0,2);
  S.Stats.RMean_std.GRF.R.Fz = [meanFz + stdFz, ...
                                meanFz, ...
                                meanFz - stdFz];

  meanMDF = mean(S.Stats.REvents.MDF);
  stdMDF = std(S.Stats.REvents.MDF);
  S.Stats.RMean_std.MDF = [meanMDF + stdMDF, meanMDF, meanMDF - stdMDF];

  meanTO = mean(S.Stats.REvents.TO);
  stdTO = std(S.Stats.REvents.TO);
  S.Stats.RMean_std.TO = [meanTO + stdTO, meanTO, meanTO - stdTO];

end
