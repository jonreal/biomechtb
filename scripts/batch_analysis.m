
% Control
if(1)
  if(~exist('c1_1'))
    c1_1 = vicon_process_trial_loc('PA_C01_SS_F0_01');
  end
  if(~exist('c1_1_mcpu'))
    c1_1_mcpu = embedded_process_data_c01('PA_C01_SS_F0_01');
  end

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

  emb_hs_timeStamp = unique(c1_1_mcpu.heelStrikeCnt)
  emb_hs_timeStamp(emb_hs_timeStamp < c1_1_mcpu.timeStamp(1)) = [];
  emb_hs_timeStamp(emb_hs_timeStamp > c1_1_mcpu.timeStamp(end)) = [];

  for i=1:numel(emb_hs_timeStamp)
    [~,ii] = min(abs(c1_1_mcpu.timeStamp - emb_hs_timeStamp(i)));
    emb_hs(i) = ii;
  end
  emb_hs(1) = [];
  emb_hs(end) = [];


  for i=1:numel(c1_1.gaitEvents.r.hs)
    [~,ii] = min(abs(c1_1_mcpu.time - c1_1.time(c1_1.gaitEvents.r.hs(i))));
    c1_1_mcpu.gaitEvents.r.hs(i) = ii;
  end

  for i=1:numel(c1_1.gaitEvents.l.hs)
    [~,ii] = min(abs(c1_1_mcpu.time - c1_1.time(c1_1.gaitEvents.l.hs(i))));
    c1_1_mcpu.gaitEvents.l.hs(i) = ii;
  end

  figure; hold all;
    plot(c1_1_mcpu.time,c1_1_mcpu.amp2s3./max(abs(c1_1_mcpu.amp2s3)),'g')
    plot(c1_1.time, ...
         c1_1.id.RAnkleMoment.X./max(abs(c1_1.id.RAnkleMoment.X)),'k')
    plot(c1_1.time(c1_1.gaitEvents.r.hs), ...
         c1_1.id.RAnkleMoment.X(c1_1.gaitEvents.r.hs) ...
          ./max(abs(c1_1.id.RAnkleMoment.X)),'or')
    plot(c1_1_mcpu.time(c1_1_mcpu.gaitEvents.r.hs),...
        c1_1_mcpu.amp2s3(c1_1_mcpu.gaitEvents.r.hs) ...
          ./max(abs(c1_1_mcpu.amp2s3)),'oy')
    title('C1-0.1 Right Heel Strike Detection','fontsize',20)

  figure; hold all;
    plot(c1_1_mcpu.time,c1_1_mcpu.amp1s3./max(abs(c1_1_mcpu.amp1s3)),'g')
    plot(c1_1.time, ...
         c1_1.id.LAnkleMoment.X./max(abs(c1_1.id.LAnkleMoment.X)),'k')
    plot(c1_1.time(c1_1.gaitEvents.l.hs), ...
         c1_1.id.LAnkleMoment.X(c1_1.gaitEvents.l.hs)...
        ./max(abs(c1_1.id.LAnkleMoment.X)),'or')
    plot(c1_1_mcpu.time(c1_1_mcpu.gaitEvents.l.hs),...
        c1_1_mcpu.amp1s3(c1_1_mcpu.gaitEvents.l.hs) ...
        ./max(abs(c1_1_mcpu.amp1s3)),'oy')
    plot(c1_1_mcpu.time(emb_hs), ...
          c1_1_mcpu.amp1s3(emb_hs)./max(abs(c1_1_mcpu.amp1s3)),'om');
    title('C1-0.1 Left Heel Strike Detection','fontsize',20)

    timeDiff = c1_1_mcpu.time(emb_hs) ...
                - c1_1_mcpu.time(c1_1_mcpu.gaitEvents.l.hs)
    avg_timeDiff = mean(timeDiff)
end

return

% Trial 0_2
if(1)
  if(~exist('t0_2'))
    t0_2 = vicon_process_trial_loc('PA_A02_SS_F0_02');
  end
  if(~exist('t0_2_mcpu'))
    t0_2_mcpu = embedded_process_data('PA_A02_SS_F0_02');
  end

  figure; hold all
    subplot(311), hold all;
    plot_std(t0_2.stats.gaitCycle,t0_2.stats.r.d_RAnkleAngles.X,0.1.*[1 1 1])
    plot_std(t0_2.stats.gaitCycle,t0_2.stats.l.d_LAnkleAngles.X,0.5*[1 1 1]);
    title({'Subject A02 - Trial-0.2', 'Velocity'},'fontsize',20)
    ylabel('rad/s','fontsize',20)
    grid on

    subplot(312), hold all;
    plot_std(t0_2.stats.gaitCycle,t0_2.stats.r.RAnkleMoment.X,0.1.*[1 1 1])
    plot_std(t0_2.stats.gaitCycle,t0_2.stats.l.LAnkleMoment.X,0.5*[1 1 1]);
    title('Moment','fontsize',20)
    ylabel('Nm/kg','fontsize',20)
    grid on

    subplot(313), hold all;
    plot_std(t0_2.stats.gaitCycle,t0_2.stats.r.RAnklePower.X,0.1.*[1 1 1])
    plot_std(t0_2.stats.gaitCycle,t0_2.stats.l.LAnklePower.X,0.5*[1 1 1]);
    title('Power','fontsize',20)
    xlabel('% Gait', 'fontsize',20)
    ylabel('W/kg','fontsize',20);
    legend('R','L')
    grid on

  for i=1:numel(t0_2.gaitEvents.r.hs)
    [~,ii] = min(abs(t0_2_mcpu.time - t0_2.time(t0_2.gaitEvents.r.hs(i))));
    t0_2_mcpu.gaitEvents.r.hs(i) = ii;
  end

  for i=1:numel(t0_2.gaitEvents.l.hs)
    [~,ii] = min(abs(t0_2_mcpu.time - t0_2.time(t0_2.gaitEvents.l.hs(i))));
    t0_2_mcpu.gaitEvents.l.hs(i) = ii;
  end

  figure; hold all;
    plot(t0_2_mcpu.time,t0_2_mcpu.amp2s3./max(abs(t0_2_mcpu.amp2s3)),'g')
    plot(t0_2.time, ...
         t0_2.id.RAnkleMoment.X./max(abs(t0_2.id.RAnkleMoment.X)),'k')
    plot(t0_2.time(t0_2.gaitEvents.r.hs), ...
         t0_2.id.RAnkleMoment.X(t0_2.gaitEvents.r.hs) ...
        ./max(abs(t0_2.id.RAnkleMoment.X)),'or')
    plot(t0_2_mcpu.time(t0_2_mcpu.gaitEvents.r.hs),...
        t0_2_mcpu.amp2s3(t0_2_mcpu.gaitEvents.r.hs)./max(abs(t0_2_mcpu.amp2s3)),'oy')
    title('A02 - 0.2 Right Heel Strike Detection','fontsize',20)

  figure; hold all;
    plot(t0_2_mcpu.time,t0_2_mcpu.amp1s3./max(abs(t0_2_mcpu.amp1s3)),'g')
    plot(t0_2.time, ...
         t0_2.id.LAnkleMoment.X./max(abs(t0_2.id.LAnkleMoment.X)),'k')
    plot(t0_2.time(t0_2.gaitEvents.l.hs), ...
         t0_2.id.LAnkleMoment.X(t0_2.gaitEvents.l.hs) ...
        ./max(abs(t0_2.id.LAnkleMoment.X)),'or')
    plot(t0_2_mcpu.time(t0_2_mcpu.gaitEvents.l.hs),...
        t0_2_mcpu.amp1s3(t0_2_mcpu.gaitEvents.l.hs)./max(abs(t0_2_mcpu.amp1s3)),'oy')
    title('A02 0.2 Left Heel Strike Detection','fontsize',20)
end

return

% Trial 1_2
if(1)
  if(~exist('t1_2'))
    t1_2 = vicon_process_trial_loc('PA_A02_SS_F1_02');
  end
  if(~exist('t1_2_mcpu'))
    t1_2_mcpu = embedded_process_data('PA_A02_SS_F1_02');
  end

  figure; hold all
    subplot(311), hold all;
    plot_std(t1_2.stats.gaitCycle,t1_2.stats.r.d_RAnkleAngles.X,0.1.*[1 1 1])
    plot_std(t1_2.stats.gaitCycle,t1_2.stats.l.d_LAnkleAngles.X,0.5*[1 1 1]);
    title({'Subject A02 - Trial-1.2', 'Velocity'},'fontsize',20)
    ylabel('rad/s','fontsize',20)
    grid on

    subplot(312), hold all;
    plot_std(t1_2.stats.gaitCycle,t1_2.stats.r.RAnkleMoment.X,0.1.*[1 1 1])
    plot_std(t1_2.stats.gaitCycle,t1_2.stats.l.LAnkleMoment.X,0.5*[1 1 1]);
    title('Moment','fontsize',20)
    ylabel('Nm/kg','fontsize',20)
    grid on

    subplot(313), hold all;
    plot_std(t1_2.stats.gaitCycle,t1_2.stats.r.RAnklePower.X,0.1.*[1 1 1])
    plot_std(t1_2.stats.gaitCycle,t1_2.stats.l.LAnklePower.X,0.5*[1 1 1]);
    title('Power','fontsize',20)
    xlabel('% Gait', 'fontsize',20)
    ylabel('W/kg','fontsize',20);
    legend('R','L')
    grid on

  emb_hs_timeStamp = unique(t1_2_mcpu.heelStrikeCnt)
  emb_hs_timeStamp(emb_hs_timeStamp < t1_2_mcpu.timeStamp(1)) = [];
  emb_hs_timeStamp(emb_hs_timeStamp > t1_2_mcpu.timeStamp(end)) = [];

  for i=1:numel(emb_hs_timeStamp)
    [~,ii] = min(abs(t1_2_mcpu.timeStamp - emb_hs_timeStamp(i)));
    emb_hs(i) = ii;
  end
  emb_hs(1) = [];
  emb_hs(end) = [];


  for i=1:numel(t1_2.gaitEvents.r.hs)
    [~,ii] = min(abs(t1_2_mcpu.time - t1_2.time(t1_2.gaitEvents.r.hs(i))));
    t1_2_mcpu.gaitEvents.r.hs(i) = ii;
  end

  for i=1:numel(t1_2.gaitEvents.l.hs)
    [~,ii] = min(abs(t1_2_mcpu.time - t1_2.time(t1_2.gaitEvents.l.hs(i))));
    t1_2_mcpu.gaitEvents.l.hs(i) = ii;
  end

  figure; hold all;
    plot(t1_2_mcpu.time,t1_2_mcpu.amp1s3./max(abs(t1_2_mcpu.amp1s3)),'g')
    plot(t1_2.time, ...
         t1_2.id.RAnkleMoment.X./max(abs(t1_2.id.RAnkleMoment.X)),'k')
    plot(t1_2.time(t1_2.gaitEvents.r.hs), ...
         t1_2.id.RAnkleMoment.X(t1_2.gaitEvents.r.hs)...
        ./max(abs(t1_2.id.RAnkleMoment.X)),'or')
    plot(t1_2_mcpu.time(t1_2_mcpu.gaitEvents.r.hs),...
        t1_2_mcpu.amp1s3(t1_2_mcpu.gaitEvents.r.hs)./max(abs(t1_2_mcpu.amp1s3)),'oy')
    title('A02 - 1.2 Right Heel Strike Detection','fontsize',20)

  figure; hold all;
    plot(t1_2_mcpu.time,t1_2_mcpu.amp2s3./max(abs(t1_2_mcpu.amp2s3)),'g')
    plot(t1_2.time, ...
         t1_2.id.LAnkleMoment.X./max(abs(t1_2.id.LAnkleMoment.X)),'k')
    plot(t1_2.time(t1_2.gaitEvents.l.hs), ...
         t1_2.id.LAnkleMoment.X(t1_2.gaitEvents.l.hs) ...
          ./max(abs(t1_2.id.LAnkleMoment.X)),'or')
    plot(t1_2_mcpu.time(t1_2_mcpu.gaitEvents.l.hs),...
        t1_2_mcpu.amp2s3(t1_2_mcpu.gaitEvents.l.hs)./max(abs(t1_2_mcpu.amp2s3)),'oy')
    plot(t1_2_mcpu.time(emb_hs), ...
          t1_2_mcpu.amp2s3(emb_hs)./max(abs(t1_2_mcpu.amp2s3)),'om');
    title('A02 1.2 Left Heel Strike Detection','fontsize',20)

    return
end

% Trial 2_2
if(1)
  if(~exist('t2_2'))
    t2_2 = vicon_process_trial_loc('PA_A02_SS_F2_02');
  end
  if(~exist('t2_2_mcpu'))
    t2_2_mcpu = embedded_process_data('PA_A02_SS_F2_02');
  end

  figure; hold all
    subplot(311), hold all;
    plot_std(t2_2.stats.gaitCycle,t2_2.stats.r.d_RAnkleAngles.X,0.1.*[1 1 1])
    plot_std(t2_2.stats.gaitCycle,t2_2.stats.l.d_LAnkleAngles.X,0.5*[1 1 1]);
    title({'Subject A02 - Trial-2.2', 'Velocity'},'fontsize',20)
    ylabel('rad/s','fontsize',20)
    grid on

    subplot(312), hold all;
    plot_std(t2_2.stats.gaitCycle,t2_2.stats.r.RAnkleMoment.X,0.1.*[1 1 1])
    plot_std(t2_2.stats.gaitCycle,t2_2.stats.l.LAnkleMoment.X,0.5*[1 1 1]);
    title('Moment','fontsize',20)
    ylabel('Nm/kg','fontsize',20)
    grid on

    subplot(313), hold all;
    plot_std(t2_2.stats.gaitCycle,t2_2.stats.r.RAnklePower.X,0.1.*[1 1 1])
    plot_std(t2_2.stats.gaitCycle,t2_2.stats.l.LAnklePower.X,0.5*[1 1 1]);
    title('Power','fontsize',20)
    xlabel('% Gait', 'fontsize',20)
    ylabel('W/kg','fontsize',20);
    legend('R','L')
    grid on

  emb_hs_timeStamp = unique(t2_2_mcpu.heelStrikeCnt)
  emb_hs_timeStamp(emb_hs_timeStamp < t2_2_mcpu.timeStamp(1)) = [];
  emb_hs_timeStamp(emb_hs_timeStamp > t2_2_mcpu.timeStamp(end)) = [];

  for i=1:numel(emb_hs_timeStamp)
    [~,ii] = min(abs(t2_2_mcpu.timeStamp - emb_hs_timeStamp(i)));
    emb_hs(i) = ii;
  end
  emb_hs(1) = [];
  emb_hs(end) = [];



  for i=1:numel(t2_2.gaitEvents.r.hs)
    [~,ii] = min(abs(t2_2_mcpu.time - t2_2.time(t2_2.gaitEvents.r.hs(i))));
    t2_2_mcpu.gaitEvents.r.hs(i) = ii;
  end

  for i=1:numel(t2_2.gaitEvents.l.hs)
    [~,ii] = min(abs(t2_2_mcpu.time - t2_2.time(t2_2.gaitEvents.l.hs(i))));
    t2_2_mcpu.gaitEvents.l.hs(i) = ii;
  end

  figure; hold all;
    plot(t2_2_mcpu.time,t2_2_mcpu.amp1s3./max(abs(t2_2_mcpu.amp1s3)),'g')
    plot(t2_2.time, ...
         t2_2.id.RAnkleMoment.X./max(abs(t2_2.id.RAnkleMoment.X)),'k')
    plot(t2_2.time(t2_2.gaitEvents.r.hs), ...
         t2_2.id.RAnkleMoment.X(t2_2.gaitEvents.r.hs) ...
        ./max(abs(t2_2.id.RAnkleMoment.X)),'or')
    plot(t2_2_mcpu.time(t2_2_mcpu.gaitEvents.r.hs),...
        t2_2_mcpu.amp1s3(t2_2_mcpu.gaitEvents.r.hs)...
        ./max(abs(t2_2_mcpu.amp1s3)),'oy')
    title('A02 - 2.2 Right Heel Strike Detection','fontsize',20)

  figure; hold all;
    plot(t2_2_mcpu.time,t2_2_mcpu.amp2s3./max(abs(t2_2_mcpu.amp2s3)),'g')
    plot(t2_2.time, ...
         t2_2.id.LAnkleMoment.X./max(abs(t2_2.id.LAnkleMoment.X)),'k')
    plot(t2_2.time(t2_2.gaitEvents.l.hs), ...
         t2_2.id.LAnkleMoment.X(t2_2.gaitEvents.l.hs) ...
        ./max(abs(t2_2.id.LAnkleMoment.X)),'or')
    plot(t2_2_mcpu.time(t2_2_mcpu.gaitEvents.l.hs),...
        t2_2_mcpu.amp2s3(t2_2_mcpu.gaitEvents.l.hs)./max(abs(t2_2_mcpu.amp2s3)),'oy')
    plot(t2_2_mcpu.time(emb_hs), ...
          t2_2_mcpu.amp2s3(emb_hs)./max(abs(t2_2_mcpu.amp2s3)),'om');
    title('A02 2.2 Left Heel Strike Detection','fontsize',20)
end

return

