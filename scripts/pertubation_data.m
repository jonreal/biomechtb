% --- Load Data ---------------------------------------------------------------
if(1)
  dataFile  = ['/Users/JonReal/Dropbox/Professional/', ...
               'UW_PHD/Prosthetic_Research/Data/PA_A01/Data/PA_C_S/', ...
               'PA_C_S.mat'];

%  speedFile  = ['/Users/JonReal/Dropbox/Professional/', ...
%               'UW_PHD/Prosthetic_Research/Data/PA_A01/Data/PA_SS/', ...
%               'PA_SS_TMf.csv'];

  load(dataFile);

%  D = csvread(speedFile,1,0);
%  tm_speed = -D(:,3);
end
% --- Analysis ----------------------------------------------------------------
time = PA_C_S.time;
RanklePos = PA_C_S.Model.RAnkleAngles.X;

COM_x = PA_C_S.Markers.CentreOfMass.X;
COM_y = PA_C_S.Markers.CentreOfMass.Y;
COM_z = PA_C_S.Markers.CentreOfMass.Z;

figure;
  subplot(411)
    plot(time,COM_x,'k')
    ylabel('$x$ (mm)','interpreter','latex','fontsize',15)
    title('Center Of Mass Trajectories', ...
          'interpreter','latex','fontsize',25)
  subplot(412)
    plot(time,COM_y,'k')
    ylabel('$y$ (mm)','interpreter','latex','fontsize',15)
  subplot(413)
    plot(time,COM_z,'k')
    ylabel('$z$ (mm)','interpreter','latex','fontsize',15)
  subplot(414)
    plot(time,PA_C_S.Model.LAnkleAngles.X,'k')
    ylabel('Angle (rad)','interpreter','latex','fontsize',15)
    xlabel('Time (s)','interpreter','latex','fontsize',15)
    ylim([-0.2, 0.5])




