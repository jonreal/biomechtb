% -- Force Sensor, gait phase detection - State Machine Script
% -
% -

dataFiles = [{'PR_SL'},{'PR_SS'}];
load_trials(dataFiles);


TS = PR_SS;
VS = PR_SL;


% -----------------------------------------------------------------------------
% Design
% -----------------------------------------------------------------------------

% Sensor Location
F1_n = 3;
F2_n = 78;
F3_n = 98;

snsr_indx = [F1_n,F2_n,F3_n];

HS = TS.Events.R.HS;
MDF = TS.Events.R.MDF;
TO = TS.Events.R.TO;

r_gaitPhase = TS.Events.R.gaitPhase;
F = [TS.Pedar.Interp.r_Pa(:,snsr_indx)];
theta = TS.Model.RAnkleAngles.X;
theta_d = TS.Model.d_RAnkleAngles.X;

%wn_Hz = 4;
%w0_Hz = 120;
%[b,a] = butter(4,wn_Hz/(w0_Hz/2));

%F_filt = filter(b,a,F)
%F_filt(F_filt < 0) = 0;

figure;
  subplot(311); hold all;
    title('Loading','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 1),'k');
    xlim([0,numel(r_gaitPhase)])
  subplot(312); hold all;
    title('Unloading','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 2)./2,'k');
    xlim([0,numel(r_gaitPhase)])
  subplot(313); hold all;
    title('Swing','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 3)./3,'k');
    xlim([0,numel(r_gaitPhase)])
    xlabel('Sample','fontsize',15)

% Heel Strike Threshold
thrshold = 25;
figure; hold all;
  title('HS Threshold','fontsize',20)
  plot(F(:,1),'k')
  plot(HS,F(HS,1),'or')
  plot(1:numel(F(:,1)),ones(size(F(:,1)))*thrshold,'--c')
  xlim([0,numel(F(:,1))]);
  grid on

% Toe Off Threshold
thrshold = 10;
figure; hold all;
  title('TO Threshold','fontsize',20)
  plot(F(:,3),'k')
  plot(TO,F(TO,3),'or')
  plot(1:numel(F(:,3)),ones(size(F(:,3)))*thrshold,'--c')
  xlim([0,numel(F(:,3))]);
  grid on

% Max Dorsiflex Threshold
thrshold = 50;
figure; hold all;
  title('MDF Threshold','fontsize',20)
  plot(F(:,2) + F(:,3))
  plot(MDF,F(MDF,2) + F(MDF,3),'or')
  plot(HS,F(HS,2) + F(HS,3),'oc')
  plot(1:numel(F(:,1)),ones(size(F(:,1)))*thrshold,'--c')
  xlim([0,numel(F(:,1))]);
  grid on

% Max Dorsiflex Threshold
thrshold = 0.2
figure; hold all;
  title('MDF Threshold','fontsize',20)
  plot(theta)
  plot(MDF,theta(MDF),'or')
  plot(HS,theta(HS),'oc')
  xlim([0,numel(F(:,1))]);
  plot(1:numel(F(:,1)),ones(size(F(:,1)))*thrshold,'--c')
  grid on

% Max Dorsiflex Threshold
figure; hold all;
  title('MDF Threshold','fontsize',20)
  plot(theta_d)
  plot(MDF,theta_d(MDF),'or')
  plot(HS,theta_d(HS),'oc')
  xlim([0,numel(F(:,1))]);
  grid on

% -----------------------------------------------------------------------------
% Test Training Set
% -----------------------------------------------------------------------------
y0 = [0,0,0];
for i=1:numel(TS.time)
  y(i,:) = gaitPhaseStateMachine(y0,F(i,1),F(i,2),F(i,3),theta_d(i));
  y0 = y(i,:);
end

figure;
  subplot(311); hold all;
    title('Taining Set -- Loading','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 1),'k');
    plot(y(:,1),'--r')
    xlim([0,numel(r_gaitPhase)])
  subplot(312); hold all;
    title('Unloading','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 2)./2,'k');
    plot(y(:,2),'--r')
    xlim([0,numel(r_gaitPhase)])
  subplot(313); hold all;
    title('Swing','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 3)./3,'k');
    plot(y(:,3),'--r')
    xlim([0,numel(r_gaitPhase)])
    xlabel('Sample','fontsize',15)

% -----------------------------------------------------------------------------
% Test Validation Set
% -----------------------------------------------------------------------------
HS = VS.Events.R.HS;
MDF = VS.Events.R.MDF;
TO = VS.Events.R.TO;

r_gaitPhase = VS.Events.R.gaitPhase;
F = [VS.Pedar.Interp.r_Pa(:,snsr_indx)];
theta = VS.Model.RAnkleAngles.X;
theta_d = VS.Model.d_RAnkleAngles.X;


% Heel Strike Threshold
thrshold = 25;
figure; hold all;
  title('HS Threshold','fontsize',20)
  plot(F(:,1),'k')
  plot(HS,F(HS,1),'or')
  plot(1:numel(F(:,1)),ones(size(F(:,1)))*thrshold,'--c')
  xlim([0,numel(F(:,1))]);
  grid on

% Toe Off Threshold
thrshold = 10;
figure; hold all;
  title('TO Threshold','fontsize',20)
  plot(F(:,3),'k')
  plot(TO,F(TO,3),'or')
  plot(1:numel(F(:,3)),ones(size(F(:,3)))*thrshold,'--c')
  xlim([0,numel(F(:,3))]);
  grid on

% Max Dorsiflex Threshold
thrshold = 50;
figure; hold all;
  title('MDF Threshold','fontsize',20)
  plot(F(:,2) + F(:,3))
  plot(MDF,F(MDF,2) + F(MDF,3),'or')
  plot(HS,F(HS,2) + F(HS,3),'oc')
  plot(1:numel(F(:,1)),ones(size(F(:,1)))*thrshold,'--c')
  xlim([0,numel(F(:,1))]);
  grid on

% Max Dorsiflex Threshold
thrshold = 0.2
figure; hold all;
  title('MDF Threshold','fontsize',20)
  plot(theta)
  plot(MDF,theta(MDF),'or')
  plot(HS,theta(HS),'oc')
  xlim([0,numel(F(:,1))]);
  plot(1:numel(F(:,1)),ones(size(F(:,1)))*thrshold,'--c')
  grid on

% Max Dorsiflex Threshold
figure; hold all;
  title('MDF Threshold','fontsize',20)
  plot(theta_d)
  plot(MDF,theta_d(MDF),'or')
  plot(HS,theta_d(HS),'oc')
  xlim([0,numel(F(:,1))]);
  grid on



y0 = [0,0,0];
for i=1:numel(VS.time)
  y(i,:) = gaitPhaseStateMachine(y0,F(i,1),F(i,2),F(i,3),theta_d(i));
  y0 = y(i,:);
end

figure;
  subplot(311); hold all;
    title('Validation Set -- Loading','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 1),'k');
    plot(y(:,1),'--r')
    xlim([0,numel(r_gaitPhase)])
  subplot(312); hold all;
    title('Unloading','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 2)./2,'k');
    plot(y(:,2),'--r')
    xlim([0,numel(r_gaitPhase)])
  subplot(313); hold all;
    title('Swing','fontsize',20)
    plot(r_gaitPhase.*(r_gaitPhase == 3)./3,'k');
    plot(y(:,3),'--r')
    xlim([0,numel(r_gaitPhase)])
    xlabel('Sample','fontsize',15)
