function load_trials(dataFile)

directory = '~/Dropbox/Professional/UW_PHD/research/data/PA_A01/data/';

  currentVars = evalin('base', 'whos');
  for i=1:numel(dataFile)
    if ismember(dataFile{i}, [currentVars(:).name])
      continue;
    end
    evalin('base', strcat('load(''', ...
          [directory,dataFile{i},'/',dataFile{i},'.mat'], ''')'));
  end
end

