  if(~exist('t0_2'))
    t0_2 = vicon_process_trial_loc('PA_A02_SS_F0_02')
  end
  if(~exist('t0_2_bbb'))
    t0_2_bbb = readData('PA_A02_SS_F0_02');
  end



  bbb_anklePos = resample(deg2rad(t0_2_bbb.anklePos./100),120,1000)
  bbb_anklePos(end) = [];
  time = [0:numel(bbb_anklePos)].* (1/120);
  time(end) = [];

  si = numel(bbb_anklePos) - numel(t0_2.id.LAnkleAngles.X) + 1;


  fobj =@(c) norm(bbb_anklePos(si:end) + c - t0_2.id.LAnkleAngles.X);
  c = fminunc(fobj,0.2178)

    figure; hold all;
      plot(t0_2.time,...
          rad2deg(t0_2.id.LAnkleAngles.X),'k')
      plot(t0_2_bbb.time ,...
           t0_2_bbb.anklePos/100 + rad2deg(c))


