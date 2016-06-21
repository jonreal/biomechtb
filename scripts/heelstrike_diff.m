if(1)

  if(~exist('c1_1'))
    c1_1 = vicon_process_trial_loc('PA_C01_SS_F0_01');
  end
  if(~exist('c1_1_mcpu'))
    c1_1_mcpu = embedded_process_data_c01('PA_C01_SS_F0_01');
  end

  fc = 60;
  fs = 1000;

  [b,a] = butter(1,2*fc/fs);
  time = c1_1_mcpu.time;
  sample = c1_1_mcpu.amp1.s3;
  [num,den] = tfdata(sys_d);
  num = num{1}
  den = den{1}

  gt = 0.*sample;
  filt = 0.*sample;
  vel = 0.*sample;

  tol = 300;
  tol2 = 450;
  % Algorithm
  for i=1:numel(time)
    if i<3
      continue
    end

    filt(i) = num(1)*sample(i) + num(2)*sample(i-1) + num(3)*sample(i-2) ...
              - den(2)*filt(i-1) - den(3)*filt(i-2);
    vel(i) = 2/Ts*filt(i) - 2/Ts*filt(i-1) - vel(i-1);

    if gt(i-1) == 0
      if (vel(i) > tol) && (filt(i) < tol2)
        gt(i) = 1;
      else
        gt(i) = 0;
      end
    elseif gt(i-1) == 1
      if (vel(i) < -tol) && (filt(i) > tol2)
        gt(i) = 0;
      else
        gt(i) = 1;
      end
    end
  end

  figure; hold all;
    plot(vel./max(abs(vel)));
    plot(sample./max(abs(sample)));
    plot(filt./max(abs(filt)));
    plot(vel.*0 + tol./max(abs(vel)));
    plot(vel.*0 - tol./max(abs(vel)));

    pause

  hs_ii = [];
  for i=2:numel(gt)
    if gt(i-1) == 0 && gt(i) == 1
      hs_ii(end+1) = i;
    end
  end
  hs_ii(1) = [];
  hs_ii(end) = [];

  figure; hold all;
    plot(sample,'k')
    plot(gt*500,'--r')
    plot(hs_ii,gt(hs_ii)*500,'om')
    grid on

  figure; hold all
    subplot(311), hold all;
    plot_std(c1_1.stats.gaitCycle,c1_1.stats.r.d_RAnkleAngles.X,0.1.*[1 1 1]);
    plot_std(c1_1.stats.gaitCycle,c1_1.stats.l.d_LAnkleAngles.X,0.5*[1 1 1]);
    title({'Subject C01 - Trial-0.1', 'Velocity'},'fontsize',20)
    ylabel('rad/s','fontsize',20)
    grid on

    subplot(312), hold all;
    plot_std(c1_1.stats.gaitCycle,c1_1.stats.r.RAnkleMoment.X,0.1.*[1 1 1])
    plot_std(c1_1.stats.gaitCycle,c1_1.stats.l.LAnkleMoment.X,0.5*[1 1 1]);
    title('Moment','fontsize',20)
    ylabel('Nm/kg','fontsize',20)
    grid on

    subplot(313), hold all;
    plot_std(c1_1.stats.gaitCycle,c1_1.stats.r.RAnklePower.X,0.1.*[1 1 1])
    plot_std(c1_1.stats.gaitCycle,c1_1.stats.l.LAnklePower.X,0.5*[1 1 1]);
    title('Power','fontsize',20)
    xlabel('% Gait', 'fontsize',20)
    ylabel('W/kg','fontsize',20);
    legend('R','L')
    grid on


  if(1)
  emb_hs_timeStamp = unique(c1_1_mcpu.heelStrikeCnt)
  emb_hs_timeStamp(emb_hs_timeStamp < c1_1_mcpu.timeStamp(1)) = [];
  emb_hs_timeStamp(emb_hs_timeStamp > c1_1_mcpu.timeStamp(end)) = [];

  for i=1:numel(emb_hs_timeStamp)
    [~,ii] = min(abs(c1_1_mcpu.timeStamp - emb_hs_timeStamp(i)));
    emb_hs(i) = ii;
  end
  emb_hs(1) = [];
  emb_hs(end) = [];

  %hs_ii = emb_hs;
  end

  for i=1:numel(c1_1.gaitEvents.r.hs)
    [~,ii] = min(abs(c1_1_mcpu.time - c1_1.time(c1_1.gaitEvents.r.hs(i))));
    c1_1_mcpu.gaitEvents.r.hs(i) = ii;
  end

  for i=1:numel(c1_1.gaitEvents.l.hs)
    [~,ii] = min(abs(c1_1_mcpu.time - c1_1.time(c1_1.gaitEvents.l.hs(i))));
    c1_1_mcpu.gaitEvents.l.hs(i) = ii;
  end

  figure; hold all;
    plot(c1_1_mcpu.time,c1_1_mcpu.amp2.s3./max(abs(c1_1_mcpu.amp2.s3)),'g')
    plot(c1_1.time, ...
         c1_1.id.RAnkleMoment.X./max(abs(c1_1.id.RAnkleMoment.X)),'k')
    plot(c1_1.time(c1_1.gaitEvents.r.hs), ...
         c1_1.id.RAnkleMoment.X(c1_1.gaitEvents.r.hs) ...
          ./max(abs(c1_1.id.RAnkleMoment.X)),'or')
    plot(c1_1_mcpu.time(c1_1_mcpu.gaitEvents.r.hs),...
        c1_1_mcpu.amp2.s3(c1_1_mcpu.gaitEvents.r.hs) ...
          ./max(abs(c1_1_mcpu.amp2.s3)),'oy')
    title('C1-0.1 Right Heel Strike Detection','fontsize',20)
    grid on;

  figure; hold all;
    plot(c1_1_mcpu.time,c1_1_mcpu.amp1.s3./max(abs(c1_1_mcpu.amp1.s3)),'g')
    plot(c1_1.time, ...
         c1_1.id.LAnkleMoment.X./max(abs(c1_1.id.LAnkleMoment.X)),'k')
    plot(c1_1.time(c1_1.gaitEvents.l.hs), ...
         c1_1.id.LAnkleMoment.X(c1_1.gaitEvents.l.hs)...
        ./max(abs(c1_1.id.LAnkleMoment.X)),'or')
    plot(c1_1_mcpu.time(c1_1_mcpu.gaitEvents.l.hs),...
        c1_1_mcpu.amp1.s3(c1_1_mcpu.gaitEvents.l.hs) ...
        ./max(abs(c1_1_mcpu.amp1.s3)),'oy')
    plot(c1_1_mcpu.time(hs_ii), ...
          c1_1_mcpu.amp1.s3(hs_ii)./max(abs(c1_1_mcpu.amp1.s3)),'om');
    title('C1-0.1 Left Heel Strike Detection','fontsize',20)
    grid on

    timeDiff = c1_1_mcpu.time(hs_ii) ...
                - c1_1_mcpu.time(c1_1_mcpu.gaitEvents.l.hs)
    avg_timeDiff = mean(timeDiff)
end


