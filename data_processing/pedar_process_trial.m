function PRSR = pedar_process_trial(trialName, sensorType)

  % Load data file
  directory = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/Data/'];
  RAW = dlmread([directory,trialName,'/',trialName,'_Pedar.asc'],'\t',7,0);

  % Load Sensor Template
  load(['~/Research_link/Prosthetic_Research/', ...
        'Source_v2/Biomechanics_Analysis_Toolbox/', ...
        'pedar_insole_ui/insole_mat_files/', ...
        sensorType,'_insole.mat']);

  % Store data
  l_x_centroid = P.centroid(1:99,1)' + 75;
  l_y_centroid = P.centroid(1:99,2)';

  r_x_centroid = P.centroid((1:99) + 100,1)' - 75;
  r_y_centroid = P.centroid((1:99) + 100,2)';

  PRSR.time = RAW(:,1) - RAW(1,1);
  PRSR.l_Pa = RAW(:,2:100);
  PRSR.r_Pa = RAW(:,101:end-1);
  PRSR.sensorTemplate.l_x_centroid = l_x_centroid;
  PRSR.sensorTemplate.l_y_centroid = l_y_centroid;
  PRSR.sensorTemplate.r_x_centroid = r_x_centroid;
  PRSR.sensorTemplate.r_y_centroid = r_y_centroid;

  PRSR.sensorTemplate.left_boundary = ...
      [P.left_boundary(:,1)+75,P.left_boundary(:,2)];
  PRSR.sensorTemplate.right_boundary = ...
      [P.right_boundary(:,1)-75,P.right_boundary(:,2)];

  figure; hold all;
    subplot(211); hold all;
      title('Right Pressure (Sensors [1,50,99])','fontsize',20);
      plot(PRSR.r_Pa(:,[1,50,99]))
   subplot(212); hold all;
      title('Left Pressure (Sensors [1,50,99])','fontsize',20);
      plot(PRSR.l_Pa(:,[1,50,99]))
end
