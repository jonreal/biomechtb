function S = bb_proces_trial(filename)

D = dlmread(filename, ' ',6,0);

S.time = D(:,1);
S.loop_time = D(:,2);
S.ankle_pos = D(:,3);
S.ankle_pos_d = D(:,4);
S.cmd = D(:,5);
S.Kp = D(:,6);
end
