function rtn = vicon_process_events(trialName,timeVector)

  % Parse model file
  file = ['./',trialName,'_Events.csv'];
  if exist(file,'file') ~= 2
    fprintf('\n\tEvent file not found: %s\n', file);
    rtn = []
    return
  end

  % Load data file
  fid = fopen(file,'r');

  % skip two words, go next line, grab sample freq
  freq = fscanf(fid, '%*s\n %f', 1);

  % skip line
  tline = fgetl(fid);

  % skip line
  tline = fgetl(fid);

  formatSpec = ['%s %s %s %f %s\n'];

  % Heel strike and toe off
  l.hs_time = [];
  l.to_time = [];
  r.hs_time = [];
  r.to_time = [];

  % Parse event file, and store event indices
  while(1)
    C = textscan(fid, formatSpec, 'delimiter', ',');
    if isempty(C{4})
      break;
    end
    if strcmp(C{2}, 'Left')
      if strcmp(C{3}, 'Foot Off')
        l.to_time(end+1,:) = C{4};
      else
        l.hs_time(end+1,:) = C{4};
      end
    else
      if strcmp(C{3}, 'Foot Off')
        r.to_time(end+1,:) = C{4};
      else
        r.hs_time(end+1,:) = C{4};
      end
    end
  end
  fclose(fid);

  % Remove nan, sort
  l.to_time = sort(l.to_time(~isnan(l.to_time)));
  l.hs_time = sort(l.hs_time(~isnan(l.hs_time)));
  r.to_time = sort(r.to_time(~isnan(r.to_time)));
  r.hs_time = sort(r.hs_time(~isnan(r.hs_time)));

  % Remove first heel strike
  l.hs_time(1) = [];
  r.hs_time(1) = [];

  % If first to is before first hs remove
  if(l.to_time(1) < l.hs_time(1))
    l.to_time(1) = [];
  end
  if(r.to_time(1) < r.hs_time(1))
    r.to_time(1) = [];
  end

  % If last to happens before last hs remove
  if(l.to_time(end) < l.hs_time(end))
    l.hs_time(end) = [];
  end
  if(r.to_time(end) < r.hs_time(end))
    r.hs_time(end) = [];
  end

  % Allocate vector size for vicon gaitevent indices
  l.hs_index = 0.*l.hs_time;
  l.to_index = 0.*l.to_time;
  r.hs_index = 0.*r.hs_time;
  r.to_index = 0.*r.to_time;

  % Convert vicon events to indices
  for i=1:numel(l.hs_time)
    [~,ii] = min(abs(timeVector - l.hs_time(i)));
    l.hs_index(i) = ii;
  end
  for i=1:numel(l.to_time)
    [~,ii] = min(abs(timeVector - l.to_time(i)));
    rtn.gaitEvents.l.to_index(i) = ii;
  end
  for i=1:numel(r.hs_time)
    [~,ii] = min(abs(timeVector - r.hs_time(i)));
    r.hs_index(i) = ii;
  end
  for i=1:numel(r.to_time)
    [~,ii] = min(abs(timeVector - r.to_time(i)));
    r.to_index(i) = ii;
  end

  gaitEvents.l = l;
  gaitEvents.r = r;

  rtn = gaitEvents;
end
