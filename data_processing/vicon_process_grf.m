function S = vicon_process_grf(trialName)

  directory = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/Data/'];
  fid = fopen([directory,trialName,'/',trialName,'_GRF.csv']);

  % skip line
  tline = fgetl(fid);

  S.freq =fscanf(fid,'%f', 1);

  S.Raw.time = [0:(1/S.freq):(30-1/S.freq)]';

  fclose(fid);

  RAW = csvread([directory,trialName,'/',trialName,'_GRF.csv'],5,0);
  S.Raw.R.Fx = RAW(:,3);
  S.Raw.R.Fy = RAW(:,4);
  S.Raw.R.Fz = RAW(:,5);

  S.Raw.L.Fx = RAW(:,6);
  S.Raw.L.Fy = RAW(:,7);
  S.Raw.L.Fz = RAW(:,8);

  % Down Sample to 120 Hz
  S.Decimate.R.Fx = decimate(RAW(:,3),10,'iir');
  S.Decimate.R.Fy = decimate(RAW(:,4),10,'iir');
  S.Decimate.R.Fz = decimate(RAW(:,5),10,'iir');

  S.Decimate.L.Fx = decimate(RAW(:,6),10,'iir');
  S.Decimate.L.Fy = decimate(RAW(:,7),10,'iir');
  S.Decimate.L.Fz = decimate(RAW(:,8),10,'iir');

  figure;
    subplot(211), hold all;
      title('GRF - Right','fontsize',20);
      plot(S.Raw.time,S.Raw.R.Fz,'o')
      stem((0:(1/120):(30-(1/120))),S.Decimate.R.Fz)
    subplot(212), hold all;
      title('GRF - Left','fontsize',20);
      plot(S.Raw.time,S.Raw.L.Fz,'o')
      stem((0:(1/120):(30-(1/120))),S.Decimate.L.Fz)
end
