% --- Script to process and save all Data trials
%
% Run after changing any of the vicon process functions.

directory = ['~/Dropbox/Professional/UW_PHD/Prosthetic_Research', ...
              '/Data/PA_A01/Data/'];

files = dir(directory);

% --- Get rid of hidden files
files = files(4:end);

for i=13:numel(files)
  close all
  fprintf(['Processing trial: ',files(i).name,'\n']);

  cmd = ['vicon_process_trial(','''',files(i).name,'''',');'];
  eval([files(i).name,' = ',cmd])

  fileName = [directory,files(i).name,'/',files(i).name,'.mat'];
  save(fileName,files(i).name);

  fprintf(['Saved ',fileName,'\n']);
  fprintf('Press Enter to continue... \n\n');
  pause
end

