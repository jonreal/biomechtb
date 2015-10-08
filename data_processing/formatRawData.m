%------------------------------------------------------------------------------
% Format Raw Data
%
% Script to convert raw gait data and format it.
%
% 19 - March - 2013
%
% Prosthetic Research v4
%------------------------------------------------------------------------------
clear all

%------------------------------------------------------------------------------
% INPUTS
%------------------------------------------------------------------------------
subjectNum = 'C14';

saveData = false;
%------------------------------------------------------------------------------
% Kinetics
%------------------------------------------------------------------------------
filename = [ 'Kinetics/',subjectNum,'.csv'];
currentLine = 1;
fid = fopen(filename);

% Crap
formatspec = '%s';
textscan(fid, formatspec, 1);
currentLine = currentLine + 1;

% Sample Frq
formatspec = '%f';
freq = textscan(fid, formatspec, 1);
currentLine = currentLine + 1;
freq = freq{:};

% Event Headers
formatspec = '%s %s %s %s %s';
eventHeaders = textscan(fid, formatspec, 1, 'delimiter', ',');
currentLine = currentLine + 1;

% Event Data
formatspec = '%s %s %s %f %s';
gaitEvents = {};
while true
  gaitEvents{end+1} = textscan(fid, formatspec, 1, 'delimiter', ',');
  currentLine = currentLine + 1;
  if (strcmp(gaitEvents{end}{5}, '')), break, end;
end

%------------------------------------------------------------------------------
% Gait Data [Angle Moment Power]
%------------------------------------------------------------------------------
% Move current line index to where data starts (4 lines after events)
currentLine = currentLine + 4;

gaitData = csvread(filename, currentLine, 0);

fclose(fid);

%------------------------------------------------------------------------------
% Find Event Indexs
%------------------------------------------------------------------------------
L_heelStrike = [];
L_toeOff = [];
R_heelStrike = [];
R_toeOff = [];
for i=1:length(gaitEvents)
  if strcmp('Foot Strike', gaitEvents{i}{3})
    if strcmp('Left', gaitEvents{i}{2})
      L_heelStrike(end+1) = floor(gaitEvents{i}{4}*freq);
    else
      R_heelStrike(end+1) = floor(gaitEvents{i}{4}*freq);
    end
  end
  if strcmp('Foot Off', gaitEvents{i}{3})
    if strcmp('Left', gaitEvents{i}{2})
      L_toeOff(end+1) = floor(gaitEvents{i}{4}*freq);
    else
      R_toeOff(end+1) = floor(gaitEvents{i}{4}*freq);
    end
  end
end
L_heelStrike = sort(L_heelStrike);
R_heelStrike = sort(R_heelStrike);
L_toeOff = sort(L_toeOff);
R_toeOff = sort(R_toeOff);

%------------------------------------------------------------------------------
% Interperolate
%------------------------------------------------------------------------------
time = gaitData(:,1)/freq;

numOfLeftCycles = length(L_heelStrike)-1;
numOfRightCycles = length(R_heelStrike)-1;

gaitCycle = linspace(0,100,1001)';
dg = gaitCycle(2) - gaitCycle(1);

L_traj = [];
for i=1:numOfLeftCycles

  % Indices
  startIndx = L_heelStrike(i);
  endIndx = L_heelStrike(i+1);

  % Raw Data
  rawTime = time(startIndx:endIndx);
  rawGaitCycle = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;i
  rawPos = deg2rad(gaitData(startIndx:endIndx,3));
  rawMom = gaitData(startIndx:endIndx,6)/1000;
  rawPow = gaitData(startIndx:endIndx,9);

  % Time Derivatives of Position
  dt = (rawTime(end) - rawTime(1))/numel(gaitCycle);
  pos = interp1(rawGaitCycle, rawPos, gaitCycle, 'spline');

  vel = gradient(pos)/dt;
  acc = gradient(vel)/dt;

  % Filter Moment, Time Derivative of Moment
  [B, A] = butter(3, 0.18);
  filtMom = filtfilt(B,A,rawMom);
  mom = interp1(rawGaitCycle, filtMom, gaitCycle, 'spline');

  ddtmom = derivative(mom)/dt;
  d2d2tmom = derivative(ddtmom)/dt;

  % Power
  pow = interp1(rawGaitCycle, rawPow, gaitCycle, 'spline');

  % Check Calculations
  if(0)
    figure;
      plot(rawGaitCycle, rawMom, rawGaitCycle, filtMom);
      title([num2str(i), ' Moment Compare']); legend('raw', 'filtered');
    figure;
      plot(gaitCycle, mom, gaitCycle, mom(1) + cumtrapz(ddtmom)*dt, '--')
      title([num2str(i), ' Moment Derivative']); legend('Mom', 'Intergral');
    figure;
      plot(gaitCycle, pos, gaitCycle, pos(1) + cumtrapz(vel)*dt, '--')
      title([num2str(i), ' Position Compare']); legend('Pos', 'Intergral');
    figure;
      plot(gaitCycle, pow, gaitCycle, -mom.*vel),
      title([num2str(i), ' Power Compare']); legend('Reported', 'Calculated');
    pause
  end

  % Pack in Matrix
  L_traj.dt(i) = dt;
  L_traj.pos(:,i) = pos;
  L_traj.vel(:,i) = vel;
  L_traj.acc(:,i) = acc;
  L_traj.mom(:,i) = mom;
  L_traj.ddtmom(:,i) = ddtmom;
  L_traj.d2d2tmom(:,i) = d2d2tmom;
  L_traj.pow(:,i) = -mom.*vel;
  L_traj.pow_report(:,i) = pow;

end

R_traj = [];
for i=1:numOfRightCycles
  % Indices
  startIndx = R_heelStrike(i);
  endIndx = R_heelStrike(i+1);

  % Raw Data
  rawTime = time(startIndx:endIndx);
  rawGaitCycle = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;
  rawPos = deg2rad(gaitData(startIndx:endIndx,12));
  rawMom = gaitData(startIndx:endIndx,15)/1000;
  rawPow = gaitData(startIndx:endIndx,18);

  % Time Derivatives of Position
  dt = (rawTime(end) - rawTime(1))/numel(gaitCycle);
  pos = interp1(rawGaitCycle, rawPos, gaitCycle, 'spline');

  vel = derivative(pos)/dt;
  acc = derivative(vel)/dt;

  % Filter Moment, Time Derivative of Moment
  [B, A] = butter(3, 0.18);
  filtMom = filtfilt(B,A,rawMom);
  mom = interp1(rawGaitCycle, filtMom, gaitCycle, 'spline');

  ddtmom = derivative(mom)/dt;
  d2d2tmom = derivative(ddtmom)/dt;

  % Power
  pow = interp1(rawGaitCycle, rawPow, gaitCycle, 'spline');

  % Check Calculations
  if(0)
    figure;
      plot(rawGaitCycle, rawMom, rawGaitCycle, filtMom);
      title([num2str(i), ' Moment Compare']); legend('raw', 'filtered');
    figure;
      plot(gaitCycle, mom, gaitCycle, mom(1) + cumtrapz(ddtmom)*dt, '--')
      title([num2str(i), ' Moment Derivative']); legend('Mom', 'Intergral');
    figure;
      plot(gaitCycle, pos, gaitCycle, pos(1) + cumtrapz(vel)*dt, '--')
      title([num2str(i), ' Position Compare']); legend('Pos', 'Intergral');
    figure;
      plot(gaitCycle, pow, gaitCycle, -mom.*vel),
      title([num2str(i), ' Power Compare']); legend('Reported', 'Calculated');
    pause
  end

  % Pack in Matrix
  R_traj.dt(i) = dt;
  R_traj.pos(:,i) = pos;
  R_traj.vel(:,i) = vel;
  R_traj.acc(:,i) = acc;
  R_traj.mom(:,i) = mom;
  R_traj.ddtmom(:,i) = ddtmom;
  R_traj.d2d2tmom(:,i) = d2d2tmom;
  R_traj.pow(:,i) = -mom.*vel;
  R_traj.pow_report(:,i) = pow;
end

%------------------------------------------------------------------------------
% Outliner for C16
%------------------------------------------------------------------------------
if strcmp(subjectNum, 'C16')
  figure;
    plot(R_traj.mom(:,10))

    R_traj.dt(10) = [];
    R_traj.pos(:,10) = [];
    R_traj.vel(:,10) = [];
    R_traj.acc(:,10) = [];
    R_traj.mom(:,10) = [];
    R_traj.ddtmom(:,10) = [];
    R_traj.d2d2tmom(:,10) = [];
    R_traj.pow(:,10) = [];
    R_traj.pow_report(:,10) = [];
end

%------------------------------------------------------------------------------
% Left Stats
%------------------------------------------------------------------------------
L_traj.dt_mean = mean(L_traj.dt);
L_traj.dt_std = std(L_traj.dt)';
L_traj.dt_pos_std = L_traj.dt_mean + L_traj.dt_std;
L_traj.dt_neg_std = L_traj.dt_mean - L_traj.dt_std;

L_traj.pos_mean = mean(L_traj.pos')';
L_traj.pos_std = std(L_traj.pos')';
L_traj.pos_pos_std = L_traj.pos_mean + L_traj.pos_std;
L_traj.pos_neg_std = L_traj.pos_mean - L_traj.pos_std;

L_traj.vel_mean = mean(L_traj.vel')';
L_traj.vel_std = std(L_traj.vel')';
L_traj.vel_pos_std = L_traj.vel_mean + L_traj.vel_std;
L_traj.vel_neg_std = L_traj.vel_mean - L_traj.vel_std;

L_traj.acc_mean = mean(L_traj.acc')';
L_traj.acc_std = std(L_traj.acc')';
L_traj.acc_pos_std = L_traj.acc_mean + L_traj.acc_std;
L_traj.acc_neg_std = L_traj.acc_mean - L_traj.acc_std;

L_traj.mom_mean = mean(L_traj.mom')';
L_traj.mom_std = std(L_traj.mom')';
L_traj.mom_pos_std = L_traj.mom_mean + L_traj.mom_std;
L_traj.mom_neg_std = L_traj.mom_mean - L_traj.mom_std;

L_traj.ddtmom_mean = mean(L_traj.ddtmom')';
L_traj.ddtmom_std = std(L_traj.ddtmom')';
L_traj.ddtmom_pos_std = L_traj.ddtmom_mean + L_traj.ddtmom_std;
L_traj.ddtmom_neg_std = L_traj.ddtmom_mean - L_traj.ddtmom_std;

L_traj.d2d2tmom_mean = mean(L_traj.d2d2tmom')';
L_traj.d2d2tmom_std = std(L_traj.d2d2tmom')';
L_traj.d2d2tmom_pos_std = L_traj.d2d2tmom_mean + L_traj.d2d2tmom_std;
L_traj.d2d2tmom_neg_std = L_traj.d2d2tmom_mean - L_traj.d2d2tmom_std;

L_traj.pow_mean = mean(L_traj.pow')';
L_traj.pow_std = std(L_traj.pow')';
L_traj.pow_pos_std = L_traj.pow_mean + L_traj.pow_std;
L_traj.pow_neg_std = L_traj.pow_mean - L_traj.pow_std;

L_traj.pow_report_mean = mean(L_traj.pow_report')';
L_traj.pow_report_std = std(L_traj.pow_report')';
L_traj.pow_report_std_pos_std = L_traj.pow_report_mean ...
                              + L_traj.pow_report_std;
L_traj.pow_report_std_neg_std = L_traj.pow_report_mean ...
                              - L_traj.pow_report_std;

% Check Calculations
if(0)
  figure;
    plot(gaitCycle, L_traj.mom_mean, ...
         gaitCycle, L_traj.mom_mean(1) ...
                    + cumtrapz(L_traj.ddtmom_mean)*L_traj.dt_mean, '--')
    title('Left Moment Derivative');
    legend('Mom Mean', 'Intergral Mean');
  figure;
    plot(gaitCycle, L_traj.pos_mean, ...
         gaitCycle, L_traj.pos_mean(1) ...
                    + cumtrapz(L_traj.vel_mean)*L_traj.dt_mean, '--')
    title('Left Position Compare');
    legend('Pos Mean', 'Intergral of Vel Mean');
  figure;
    plot(gaitCycle, L_traj.pow_report_mean, ...
         gaitCycle, L_traj.pow_mean),
    title('Left Power Compare');
    legend('Reported Mean', 'Calculated Mean');
end

%------------------------------------------------------------------------------
% Right Stats
%------------------------------------------------------------------------------
R_traj.dt_mean = mean(R_traj.dt);
R_traj.dt_std = std(R_traj.dt)';
R_traj.dt_pos_std = R_traj.dt_mean + R_traj.dt_std;
R_traj.dt_neg_std = R_traj.dt_mean - R_traj.dt_std;

R_traj.pos_mean = mean(R_traj.pos')';
R_traj.pos_std = std(R_traj.pos')';
R_traj.pos_pos_std = R_traj.pos_mean + R_traj.pos_std;
R_traj.pos_neg_std = R_traj.pos_mean - R_traj.pos_std;

R_traj.vel_mean = mean(R_traj.vel')';
R_traj.vel_std = std(R_traj.vel')';
R_traj.vel_pos_std = R_traj.vel_mean + R_traj.vel_std;
R_traj.vel_neg_std = R_traj.vel_mean - R_traj.vel_std;

R_traj.acc_mean = mean(R_traj.acc')';
R_traj.acc_std = std(R_traj.acc')';
R_traj.acc_pos_std = R_traj.acc_mean + R_traj.acc_std;
R_traj.acc_neg_std = R_traj.acc_mean - R_traj.acc_std;

R_traj.mom_mean = mean(R_traj.mom')';
R_traj.mom_std = std(R_traj.mom')';
R_traj.mom_pos_std = R_traj.mom_mean + R_traj.mom_std;
R_traj.mom_neg_std = R_traj.mom_mean - R_traj.mom_std;

R_traj.ddtmom_mean = mean(R_traj.ddtmom')';
R_traj.ddtmom_std = std(R_traj.ddtmom')';
R_traj.ddtmom_pos_std = R_traj.ddtmom_mean + R_traj.ddtmom_std;
R_traj.ddtmom_neg_std = R_traj.ddtmom_mean - R_traj.ddtmom_std;

R_traj.d2d2tmom_mean = mean(R_traj.d2d2tmom')';
R_traj.d2d2tmom_std = std(R_traj.d2d2tmom')';
R_traj.d2d2tmom_pos_std = R_traj.d2d2tmom_mean + R_traj.d2d2tmom_std;
R_traj.d2d2tmom_neg_std = R_traj.d2d2tmom_mean - R_traj.d2d2tmom_std;

R_traj.pow_mean = mean(R_traj.pow')';
R_traj.pow_std = std(R_traj.pow')';
R_traj.pow_pos_std = R_traj.pow_mean + R_traj.pow_std;
R_traj.pow_neg_std = R_traj.pow_mean - R_traj.pow_std;

R_traj.pow_report_mean = mean(R_traj.pow_report')';
R_traj.pow_report_std = std(R_traj.pow_report')';
R_traj.pow_report_std_pos_std = R_traj.pow_report_mean ...
                              + R_traj.pow_report_std;
R_traj.pow_report_std_neg_std = R_traj.pow_report_mean ...
                              - R_traj.pow_report_std;

% Check Calculations
if(0)
  figure;
    plot(gaitCycle, R_traj.mom_mean, ...
         gaitCycle, R_traj.mom_mean(1) ...
                    + cumtrapz(R_traj.ddtmom_mean)*R_traj.dt_mean, '--')
    title('Right Moment Derivative');
    legend('Mom Mean', 'Intergral Mean');
  figure;
    plot(gaitCycle, R_traj.pos_mean, ...
         gaitCycle, R_traj.pos_mean(1) ...
                    + cumtrapz(R_traj.vel_mean)*R_traj.dt_mean, '--')
    title('Right Position Compare');
    legend('Pos Mean', 'Intergral of Vel Mean');
  figure;
    plot(gaitCycle, R_traj.pow_report_mean, ...
         gaitCycle, R_traj.pow_mean),
    title('Right Power Compare');
    legend('Reported Mean', 'Calculated Mean');
end

%------------------------------------------------------------------------------
% Combined Stats
%------------------------------------------------------------------------------
Both_traj.dt_mean = mean([L_traj.dt, R_traj.dt]')';
Both_traj.dt_std = std([L_traj.dt R_traj.dt]')';
Both_traj.dt_pos_std = Both_traj.dt_mean + Both_traj.dt_std;
Both_traj.dt_neg_std = Both_traj.dt_mean - Both_traj.dt_std;

Both_traj.pos_mean = mean([L_traj.pos, R_traj.pos]')';
Both_traj.pos_std = std([L_traj.pos, R_traj.pos]')';
Both_traj.pos_pos_std = Both_traj.pos_mean + Both_traj.pos_std;
Both_traj.pos_neg_std = Both_traj.pos_mean - Both_traj.pos_std;

Both_traj.vel_mean = mean([L_traj.vel, R_traj.vel]')';
Both_traj.vel_std = std([L_traj.vel, R_traj.vel]')';
Both_traj.vel_pos_std = Both_traj.vel_mean + Both_traj.vel_std;
Both_traj.vel_neg_std = Both_traj.vel_mean - Both_traj.vel_std;

Both_traj.acc_mean = mean([L_traj.acc, R_traj.acc]')';
Both_traj.acc_std = std([L_traj.acc, R_traj.acc]')';
Both_traj.acc_pos_std = Both_traj.acc_mean + Both_traj.acc_std;
Both_traj.acc_neg_std = Both_traj.acc_mean - Both_traj.acc_std;

Both_traj.mom_mean = mean([L_traj.mom, R_traj.mom]')';
Both_traj.mom_std = std([L_traj.mom, R_traj.mom]')';
Both_traj.mom_pos_std = Both_traj.mom_mean + Both_traj.mom_std;
Both_traj.mom_neg_std = Both_traj.mom_mean - Both_traj.mom_std;

Both_traj.ddtmom_mean = mean([L_traj.ddtmom, R_traj.ddtmom]')';
Both_traj.ddtmom_std = std([L_traj.ddtmom, R_traj.ddtmom]')';
Both_traj.ddtmom_pos_std = Both_traj.ddtmom_mean + Both_traj.ddtmom_std;
Both_traj.ddtmom_neg_std = Both_traj.ddtmom_mean - Both_traj.ddtmom_std;

Both_traj.d2d2tmom_mean = mean([L_traj.d2d2tmom, R_traj.d2d2tmom]')';
Both_traj.d2d2tmom_std = std([L_traj.d2d2tmom, R_traj.d2d2tmom]')';
Both_traj.d2d2tmom_pos_std = Both_traj.d2d2tmom_mean + Both_traj.d2d2tmom_std;
Both_traj.d2d2tmom_neg_std = Both_traj.d2d2tmom_mean - Both_traj.d2d2tmom_std;

Both_traj.pow_mean = mean([L_traj.pow, R_traj.pow]')';
Both_traj.pow_std = std([L_traj.pow, R_traj.pow]')';
Both_traj.pow_pos_std = Both_traj.pow_mean + Both_traj.pow_std;
Both_traj.pow_neg_std = Both_traj.pow_mean - Both_traj.pow_std;

Both_traj.pow_report_mean = mean([L_traj.pow_report, R_traj.pow_report]')';
Both_traj.pow_report_std = std([L_traj.pow_report, R_traj.pow_report]')';
Both_traj.pow_report_std_pos_std = Both_traj.pow_report_mean ...
                              + Both_traj.pow_report_std;
Both_traj.pow_report_std_neg_std = Both_traj.pow_report_mean ...
                              - Both_traj.pow_report_std;

% Check Calculations
if(0)
  figure;
    plot(gaitCycle, Both_traj.mom_mean, ...
         gaitCycle, Both_traj.mom_mean(1) ...
                    + cumtrapz(Both_traj.ddtmom_mean)*Both_traj.dt_mean, '--')
    title('Combined Moment Derivative');
    legend('Mom Mean', 'Intergral Mean');
  figure;
    plot(gaitCycle, Both_traj.pos_mean, ...
         gaitCycle, Both_traj.pos_mean(1) ...
                    + cumtrapz(Both_traj.vel_mean)*Both_traj.dt_mean, '--')
    title('Combined Position Compare');
    legend('Pos Mean', 'Intergral of Vel Mean');
  figure;
    plot(gaitCycle, Both_traj.pow_report_mean, ...
         gaitCycle, Both_traj.pow_mean),
    title('Combined Power Compare');
    legend('Reported Mean', 'Calculated Mean');
end


%------------------------------------------------------------------------------
% Plot Results
%------------------------------------------------------------------------------
figure;
  subplot(321)
    plot(gaitCycle, L_traj.pos)
    title('Left')
    ylabel('Angle [rad]')

  subplot(322)
    plot(gaitCycle, R_traj.pos);
    title('Right')

  subplot(323)
    plot(gaitCycle, L_traj.mom);
    ylabel('Moment [N/kg]')

  subplot(324)
    plot(gaitCycle, R_traj.mom);

  subplot(325)
    plot(gaitCycle, L_traj.pow);
    xlabel('% gait')
    ylabel('Power [W]')

  subplot(326)
    plot(gaitCycle, R_traj.pow);
    xlabel('% gait')

figure;
  subplot(321)
    plot(gaitCycle, [L_traj.pos_neg_std, L_traj.pos_mean, L_traj.pos_pos_std]);
    title('Left')
    ylabel('Angle [rad]')

  subplot(322)
    plot(gaitCycle, [R_traj.pos_neg_std, R_traj.pos_mean, R_traj.pos_pos_std]);
    title('Right')

  subplot(323)
    plot(gaitCycle, [L_traj.mom_neg_std, L_traj.mom_mean, L_traj.mom_pos_std]);
    ylabel('Moment [N/kg]')

  subplot(324)
    plot(gaitCycle, [R_traj.mom_neg_std, R_traj.mom_mean, R_traj.mom_pos_std]);

  subplot(325)
    plot(gaitCycle, [L_traj.pow_neg_std, L_traj.pow_mean, L_traj.pow_pos_std]);
    xlabel('% gait')
    ylabel('Power [W]')

  subplot(326)
    plot(gaitCycle, [R_traj.pow_neg_std, R_traj.pow_mean, R_traj.pow_pos_std]);
    xlabel('% gait')


figure;
  subplot(2,1,1)
    plot(-L_traj.pos_mean, L_traj.mom_mean);
  subplot(2,1,2)
    plot(-R_traj.pos_mean, R_traj.mom_mean);

figure;
  plot(-Both_traj.pos_mean, Both_traj.mom_mean);

%------------------------------------------------------------------------------
% Pack Results
%------------------------------------------------------------------------------
% Height, Weight, Speed, Age info
controlInfo = { {'C12'; 'C13'; 'C14';
                 'C15'; 'C16'; 'C18';
                 'C19'; 'C21'; 'C22';
                 'C23'; 'C24'; 'C25';},
                 [1.84, 94.1, 1.16, 30;
                  1.7, 57.3, 1.16, 28;
                  1.84, 84.1, 1.1, 57;
                  1.69, 62.7, 1.16, 23;
                  1.86, 92.3, 1.2, 59;
                  1.77, 82.3, 1.28, 29;
                  1.78, 95.9, 1.24, 29;
                  1.81, 94.5, 1.26, 42;
                  1.74, 79.1, 1.44, 48;
                  1.62, 59.5, 1.34, 35;
                  1.62, 61.8, 0.98, 63;
                  1.79, 83.6, 1.42, 57;]};


index = find(strcmp(subjectNum, controlInfo{1}) == 1);
S.height = controlInfo{2}(index,1);
S.weight = controlInfo{2}(index,2);
S.speed = controlInfo{2}(index,3);
S.age = controlInfo{2}(index,4);

S.L_traj.dt = L_traj.dt_mean;
S.L_traj.gaitCycle = gaitCycle;
S.L_traj.xd = L_traj.pos_mean;
S.L_traj.x1d = L_traj.vel_mean;
S.L_traj.x2d = L_traj.acc_mean;
S.L_traj.Td = -L_traj.mom_mean;
S.L_traj.T1d = -L_traj.ddtmom_mean;
S.L_traj.T2d = -L_traj.d2d2tmom_mean;
S.L_traj.Pd = L_traj.pow_mean;
S.L_traj.Pd_report = L_traj.pow_report_mean;

S.R_traj.dt = R_traj.dt_mean;
S.R_traj.gaitCycle = gaitCycle;
S.R_traj.xd = R_traj.pos_mean;
S.R_traj.x1d = R_traj.vel_mean;
S.R_traj.x2d = R_traj.acc_mean;
S.R_traj.Td = -R_traj.mom_mean;
S.R_traj.T1d = -R_traj.ddtmom_mean;
S.R_traj.T2d = -R_traj.d2d2tmom_mean;
S.R_traj.Pd = R_traj.pow_mean;
S.R_traj.Pd_report = R_traj.pow_report_mean;

S.Both_traj.dt = Both_traj.dt_mean;
S.Both_traj.gaitCycle = gaitCycle;
S.Both_traj.xd = Both_traj.pos_mean;
S.Both_traj.x1d = Both_traj.vel_mean;
S.Both_traj.x2d = Both_traj.acc_mean;
S.Both_traj.Td = -Both_traj.mom_mean;
S.Both_traj.T1d = -Both_traj.ddtmom_mean;
S.Both_traj.T2d = -Both_traj.d2d2tmom_mean;
S.Both_traj.Pd = Both_traj.pow_mean;
S.Both_traj.Pd_report = Both_traj.pow_report_mean;

if saveData, save(['../dataComp/', subjectNum, 'Data'], 'S'), end
