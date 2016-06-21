function S = vicon_process_model(trialName)
% ---
% --- Function processes vicon model data (Joint kinetics/kinematics and
% markers)

  debug = 0;

  directory = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/Data/'];
  fid = fopen([directory,trialName,'/',trialName,'_Model.csv']);

  % Store trialName
  S.name = trialName;

  % Hard code mass for now
  S.mass = 77.7;

  % skip two words, go next line, grab sample freq
  S.freq = fscanf(fid, '%*s %*s\n %f', 1);

  % skip line
  tline = fgetl(fid);

  % parse labels and create struct fields names
  tline = fgetl(fid);
  labels = textscan(tline, '%s', 'delimiter', ',');
  [m, n] = size(labels{1});
  labelCnt = 0;
  for i=1:m
    if ~(strcmp(labels{1}{i},''))
      % pasre header on label
      [token, remain] = strtok(labels{1}{i},':');
      labelCnt = labelCnt + 1;
      labelName{labelCnt} = remain(2:end);
    end
  end

  % find number of subfields
  tline = fgetl(fid);
  sublabels = textscan(tline, '%s', 'delimiter', ',');
  [m, n] = size(sublabels{1});

  % skip line of units
  fgetl(fid);

  % put data in matrix
  D = zeros(1,m);
  formatSpec = [repmat(['%f'],1,m-1), '%f\n'];
  rowCount = 1;
  while(1)
    C = textscan(fid, formatSpec, 'delimiter', ',');
    dataRow = cell2mat(C);
    if isempty(dataRow)
      break;
    end
    D(rowCount,:) = cell2mat(C);
    rowCount = rowCount + 1;
  end

  % make time vector
  S.time = [0:(1/S.freq):(numel(D(:,1))*(1/S.freq) - (1/S.freq))]';

  % Butterworth 4th cutoff = 6 HZ;
  [b,a] = butter(4,2*(6/120));


  % put data (X,Y,Z subfields per field)
  % Fix units
  % Pattern is (L/R)(joint)Angle, Moment, Power
  for i=1:(labelCnt/3)

    % Angle -------------------------------------------------------------------
    jointAngle_X = deg2rad(D(:,3 + (9*i - 9)));
    jointAngle_Y = deg2rad(D(:,4 + (9*i - 9)));
    jointAngle_Z = deg2rad(D(:,5 + (9*i - 9)));

    % Set NaN to zero
    jointAngle_X(isnan(jointAngle_X)) = 0;
    jointAngle_Y(isnan(jointAngle_Y)) = 0;
    jointAngle_Z(isnan(jointAngle_Z)) = 0;

    % Filter
    jointAngle_X = filtfilt(b,a,jointAngle_X);
    jointAngle_Y = filtfilt(b,a,jointAngle_Y);
    jointAngle_Z = filtfilt(b,a,jointAngle_Z);

    d_jointAngle_X = gradient(jointAngle_X)*S.freq;
    d_jointAngle_Y = gradient(jointAngle_Y)*S.freq;
    d_jointAngle_Z = gradient(jointAngle_Z)*S.freq;

    d2_jointAngle_X = gradient(d_jointAngle_X)*S.freq;
    d2_jointAngle_Y = gradient(d_jointAngle_Y)*S.freq;
    d2_jointAngle_Z = gradient(d_jointAngle_Z)*S.freq;

    % Moment ------------------------------------------------------------------
    jointMoment_X = -D(:,6 + (9*i - 9))/1000;
    jointMoment_Y = -D(:,7 + (9*i - 9))/1000;
    jointMoment_Z = -D(:,8 + (9*i - 9))/1000;

    % Set NaN to zero
    jointMoment_X(isnan(jointMoment_X))= 0;
    jointMoment_Y(isnan(jointMoment_Y))= 0;
    jointMoment_Z(isnan(jointMoment_Z))= 0;

    % Filter
    jointMoment_X = filtfilt(b,a,jointMoment_X);
    jointMoment_Y = filtfilt(b,a,jointMoment_Y);
    jointMoment_Z = filtfilt(b,a,jointMoment_Z);

    d_jointMoment_X = gradient(jointMoment_X)*S.freq;
    d_jointMoment_Y = gradient(jointMoment_Y)*S.freq;
    d_jointMoment_Z = gradient(jointMoment_Z)*S.freq;

    d2_jointMoment_X = gradient(d_jointMoment_X)*S.freq;
    d2_jointMoment_Y = gradient(d_jointMoment_Y)*S.freq;
    d2_jointMoment_Z = gradient(d_jointMoment_Z)*S.freq;

    % Power (calculate) -------------------------------------------------------
    jointPower_X = d_jointAngle_X.*jointMoment_X;
    jointPower_Y = d_jointAngle_Y.*jointMoment_Y;
    jointPower_Z = d_jointAngle_Z.*jointMoment_Z;

    d_jointPower_X = gradient(jointPower_X)*S.freq;
    d_jointPower_Y = gradient(jointPower_Y)*S.freq;
    d_jointPower_Z = gradient(jointPower_Z)*S.freq;

    d2_jointPower_X = gradient(d_jointPower_X)*S.freq;
    d2_jointPower_Y = gradient(d_jointPower_Y)*S.freq;
    d2_jointPower_Z = gradient(d_jointPower_Z)*S.freq;

    % Store Angle info --------------------------------------------------------
    S.Model.(labelName{1 + (3*i - 3)}).X = jointAngle_X;
    S.Model.(labelName{1 + (3*i - 3)}).Y = jointAngle_Y;
    S.Model.(labelName{1 + (3*i - 3)}).Z = jointAngle_Z;

    % 1st time derivative
    S.Model.(['d_',labelName{1 + (3*i - 3)}]).X = d_jointAngle_X;
    S.Model.(['d_',labelName{1 + (3*i - 3)}]).Y = d_jointAngle_Y;
    S.Model.(['d_',labelName{1 + (3*i - 3)}]).Z = d_jointAngle_Z;

    % 2nd time derivative
    S.Model.(['d2_',labelName{1 + (3*i - 3)}]).X = d2_jointAngle_X;
    S.Model.(['d2_',labelName{1 + (3*i - 3)}]).Y = d2_jointAngle_Y;
    S.Model.(['d2_',labelName{1 + (3*i - 3)}]).Z = d2_jointAngle_Z;

    if(debug)
      figure;
        subplot(311); hold all;
          title(['Saggital Plane ', labelName{1 + (3*i - 3)}],'fontsize',20)
          plot(jointAngle_X)
        subplot(312); hold all;
          plot(d_jointAngle_X)
        subplot(313); hold all;
          plot(d2_jointAngle_X)
    end

    % Store Moment info --------------------------------------------------------
    S.Model.(labelName{2 + (3*i - 3)}).X = jointMoment_X;
    S.Model.(labelName{2 + (3*i - 3)}).Y = jointMoment_Y;
    S.Model.(labelName{2 + (3*i - 3)}).Z = jointMoment_Z;

    % 1st time derivative
    S.Model.(['d_',labelName{2 + (3*i - 3)}]).X = d_jointMoment_X;
    S.Model.(['d_',labelName{2 + (3*i - 3)}]).Y = d_jointMoment_Y;
    S.Model.(['d_',labelName{2 + (3*i - 3)}]).Z = d_jointMoment_Z;

    % 2nd time derivative
    S.Model.(['d2_',labelName{2 + (3*i - 3)}]).X = d2_jointMoment_X;
    S.Model.(['d2_',labelName{2 + (3*i - 3)}]).Y = d2_jointMoment_Y;
    S.Model.(['d2_',labelName{2 + (3*i - 3)}]).Z = d2_jointMoment_Z;

    if(debug)
      figure;
        subplot(311); hold all;
          title(['Saggital Plane ', labelName{2 + (3*i - 3)}],'fontsize',20)
          plot(jointMoment_X)
        subplot(312); hold all;
          plot(d_jointMoment_X)
        subplot(313); hold all;
          plot(d2_jointMoment_X)
    end

    % Store Power info --------------------------------------------------------
    S.Model.(labelName{3 + (3*i - 3)}).X = jointPower_X;
    S.Model.(labelName{3 + (3*i - 3)}).Y = jointPower_Y;
    S.Model.(labelName{3 + (3*i - 3)}).Z = jointPower_Z;

    % 1st time derivative
    S.Model.(['d_',labelName{3 + (3*i - 3)}]).X = d_jointPower_X;
    S.Model.(['d_',labelName{3 + (3*i - 3)}]).Y = d_jointPower_Y;
    S.Model.(['d_',labelName{3 + (3*i - 3)}]).Z = d_jointPower_Z;

    % 2nd time derivative
    S.Model.(['d2_',labelName{3 + (3*i - 3)}]).X = d2_jointPower_X;
    S.Model.(['d2_',labelName{3 + (3*i - 3)}]).Y = d2_jointPower_Y;
    S.Model.(['d2_',labelName{3 + (3*i - 3)}]).Z = d2_jointPower_Z;

    if(debug)
      figure;
        subplot(311); hold all;
          title(['Saggital Plane ', labelName{3 + (3*i - 3)}],'fontsize',20)
          plot(jointPower_X)
        subplot(312); hold all;
          plot(d_jointPower_X)
        subplot(313); hold all;
          plot(d2_jointPower_X)
    end

  end

  % Marker Trajectories
  % skip line
  tline = fgetl(fid);

  % skip line
  tline = fgetl(fid);

  % parse labels and create struct fields names
  tline = fgetl(fid);
  labels = textscan(tline, '%s', 'delimiter', ',');
  [m, n] = size(labels{1});
  labelCnt = 0;
  for i=1:m
    if ~(strcmp(labels{1}{i},''))
      % pasre header on label
      [token, remain] = strtok(labels{1}{i},':');
      labelCnt = labelCnt + 1;
      labelName{labelCnt} = remain(2:end);
    end
  end

  % find number of subfields
  tline = fgetl(fid);
  sublabels = textscan(tline, '%s', 'delimiter', ',');
  [m, n] = size(sublabels{1});

  % skip line of units
  fgetl(fid);

  % put data in matrix
  D = zeros(1,m);
  formatSpec = [repmat(['%f'],1,m-1), '%f\n'];
  rowCount = 1;
  while(1)
    C = textscan(fid, formatSpec, 'delimiter', ',');
    dataRow = cell2mat(C);
    if isempty(dataRow)
      break;
    end
    D(rowCount,:) = cell2mat(C);
    rowCount = rowCount + 1;
  end

  % put data (X,Y,Z subfields per field)
  for i=1:labelCnt
    S.Markers.(labelName{i}).X = D(:,3 + (3*i - 3));
    S.Markers.(labelName{i}).Y = D(:,4 + (3*i - 3));
    S.Markers.(labelName{i}).Z = D(:,5 + (3*i - 3));
  end

  fclose(fid);
end
