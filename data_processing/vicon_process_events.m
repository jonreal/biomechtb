function Events = vicon_process_events(S)


  % Load data file
  directory = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/Data/'];
  fid = fopen([directory,S.name,'/',S.name,'_Events.csv']);

  % skip two words, go next line, grab sample freq
  freq = fscanf(fid, '%*s\n %f', 1);

  % skip line
  tline = fgetl(fid);

  % skip line
  tline = fgetl(fid);

  formatSpec = ['%s %s %s %f %s\n'];
  LFootOff_indx = [];
  LHeelStrike_indx = [];
  RFootOff_indx = [];
  RHeelStrike_indx = [];

  while(1)
    C = textscan(fid, formatSpec, 'delimiter', ',');
    if isempty(C{4})
      break;
    end
    if strcmp(C{2}, 'Left')
      if strcmp(C{3}, 'Foot Off')
        LFootOff_indx(end+1,:) = floor(C{4}*freq);
      else
        LHeelStrike_indx(end+1,:) = floor(C{4}*freq);
      end
    else
      if strcmp(C{3}, 'Foot Off')
        RFootOff_indx(end+1,:) = floor(C{4}*freq);
      else
        RHeelStrike_indx(end+1,:) = floor(C{4}*freq);
      end
    end
  end

  % Remove nan
  Events.L.TO = sort(LFootOff_indx(~isnan(LFootOff_indx)));
  Events.L.HS = sort(LHeelStrike_indx(~isnan(LHeelStrike_indx)));
  Events.R.TO = sort(RFootOff_indx(~isnan(RFootOff_indx)));
  Events.R.HS = sort(RHeelStrike_indx(~isnan(RHeelStrike_indx)));

  GRF = {S.GRF.Decimate.R.Fz,S.GRF.Decimate.L.Fz};
  POS = {S.Model.RAnkleAngles.X,S.Model.LAnkleAngles.X};
  VEL = {S.Model.d_RAnkleAngles.X,S.Model.d_LAnkleAngles.X};
  MOM = {S.Model.RAnkleMoment.X,S.Model.LAnkleMoment.X};

  FOOT = {'R','L'};

  for k=1:2

    HS = Events.(FOOT{k}).HS + 1;
    TO = Events.(FOOT{k}).TO + 1;
    HS(HS == 1) = [];
    TO(TO == 1) = [];

    % Find MDF (Peaks)
    pos = POS{k};
    vel = VEL{k};
    mom = MOM{k};

    [pks,MDF] = findpeaks(-mom);
    MDF(pks < 0.8) = [];

    figure; hold all;
      plot(mom)
      plot(MDF,mom(MDF),'or')

    % Remove false MDF's
%    ii = find(MDF(2:end) - MDF(1:end-1) < 100);
%    ii = ii(1:2:numel(ii));
%    MDF(ii+1) = [];

    % Create ground truth vector
    labels = [repmat('H',size(HS)); ...
             repmat('M',size(MDF)); ...
             repmat('T',size(TO))];
    indices = [HS; MDF; TO];

    [~,ii] = sort(indices);
    indices = indices(ii);
    labels = labels(ii);

    gaitPhase = nan*ones(size(S.time));
    if strcmp(labels(1),'H')
      gaitPhase(1:indices(1)) = 3;
      phase = 3;
    elseif strcmp(labels(1),'M')
      gaitPhase(1:indices(1)) = 1;
      phase = 1;
    elseif strcmp(labels(1),'T')
      gaitPhase(1:indices(1)) = 2;
      phase = 2;
    end

    for i=1:(numel(indices)-1)
      phase = phase + 1;
      if phase == 4;
        phase = 1;
      end
      gaitPhase(indices(i)+1:indices(i+1)) = phase;
    end

    if strcmp(labels(end),'H')
      gaitPhase(indices(end)+1:end) = 1;
    elseif strcmp(labels(end),'M')
      gaitPhase(indices(end)+1:end) = 2;
    elseif strcmp(labels(end),'T')
      gaitPhase(indices(end)+1:end) = 3;
    end


    Events.(FOOT{k}).HS = HS;
    Events.(FOOT{k}).TO = TO;
    Events.(FOOT{k}).MDF = MDF;
    Events.(FOOT{k}).gaitPhase = gaitPhase;

   end

    figure;
      subplot(211); hold all;
        plot(S.Model.RAnkleAngles.X,'k')
        [AX,H1,H2] = plotyy(1:numel(S.Model.RAnkleAngles.X), ...
                            S.Model.RAnkleAngles.X, ...
                            1:numel(Events.R.gaitPhase),...
                            Events.R.gaitPhase);
        delete(H1);
        set(H2,'Marker','o','LineStyle','none')
        plot(Events.R.HS,S.Model.RAnkleAngles.X(Events.R.HS),'or')
        plot(Events.R.TO,S.Model.RAnkleAngles.X(Events.R.TO),'og')
        plot(Events.R.MDF,S.Model.RAnkleAngles.X(Events.R.MDF),'oc')


      subplot(212); hold all;
        plot(S.Model.LAnkleAngles.X,'k')
        [AX,H1,H2] = plotyy(1:numel(S.Model.LAnkleAngles.X), ...
                            S.Model.LAnkleAngles.X, ...
                            1:numel(Events.L.gaitPhase),...
                            Events.L.gaitPhase);
        delete(H1);
        set(H2,'Marker','o','LineStyle','none')
        plot(Events.L.HS,S.Model.LAnkleAngles.X(Events.L.HS),'or')
        plot(Events.L.TO,S.Model.LAnkleAngles.X(Events.L.TO),'og')
        plot(Events.L.MDF,S.Model.LAnkleAngles.X(Events.L.MDF),'oc')
        legend('pos','HS','TO','MDF')

    figure;
      subplot(211); hold all;
        plot(S.GRF.Decimate.R.Fz,'k')
        plot(Events.R.HS,S.GRF.Decimate.R.Fz(Events.R.HS),'or')
        plot(Events.R.TO,S.GRF.Decimate.R.Fz(Events.R.TO),'og')
        plot(Events.R.MDF,S.GRF.Decimate.R.Fz(Events.R.MDF),'oc')
        title('Right GRF')

      subplot(212); hold all;
        plot(S.GRF.Decimate.L.Fz,'k')
        plot(Events.L.HS,S.GRF.Decimate.L.Fz(Events.L.HS),'or')
        plot(Events.L.TO,S.GRF.Decimate.L.Fz(Events.L.TO),'og')
        plot(Events.L.MDF,S.GRF.Decimate.L.Fz(Events.L.MDF),'oc')
        title('Left GRF')
        legend('Fz','HS','TO','MDF')
end
