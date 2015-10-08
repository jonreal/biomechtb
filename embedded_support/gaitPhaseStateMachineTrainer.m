function S = gaitPhaseStateMachineTrainer(TS)


  F1_n = 3;
  F2_n = 58;
  F3_n = 98;

  snsr_indx = [F1_n,F2_n,F3_n];

  % --- Ground Truth Transitions
  HS = TS.Events.R.HS;
  MDF = TS.Events.R.MDF;
  TO = TS.Events.R.TO;

  % --- Time Signals
  r_gaitPhase = TS.Events.R.gaitPhase;
  F1 = TS.Pedar.Interp.r_Pa(:,F1_n);
  F2 = TS.Pedar.Interp.r_Pa(:,F2_n);
  F3 = TS.Pedar.Interp.r_Pa(:,F3_n);
  theta_d = TS.Model.d_RAnkleAngles.X;

  % --- Mean Signals
  F1_mean_std = [TS.Stats.RMean_std.r_Pa{1}(:,F1_n),...
                 TS.Stats.RMean_std.r_Pa{2}(:,F1_n), ...
                 TS.Stats.RMean_std.r_Pa{3}(:,F1_n)];
  F2_mean_std = [TS.Stats.RMean_std.r_Pa{1}(:,F2_n),...
                 TS.Stats.RMean_std.r_Pa{2}(:,F2_n), ...
                 TS.Stats.RMean_std.r_Pa{3}(:,F2_n)];
  F3_mean_std = [TS.Stats.RMean_std.r_Pa{1}(:,F3_n),...
                 TS.Stats.RMean_std.r_Pa{2}(:,F3_n), ...
                 TS.Stats.RMean_std.r_Pa{3}(:,F3_n)];

  GRF_mean_std = TS.Stats.RMean_std.GRF.R.Fz;
  theta_mean_std = TS.Stats.RMean_std.RAnkleAngles.X;
  theta_d_mean_std = TS.Stats.RMean_std.d_RAnkleAngles.X;


  gaitCycle = linspace(0,100,1001);
  figure;
    subplot(611); hold all;
      plot_std(1,gaitCycle,F1_mean_std)
      plot(gaitCycle,F1_mean_std(:,2),'k')
    subplot(612); hold all;
      plot_std(1,gaitCycle,F2_mean_std)
      plot(gaitCycle,F2_mean_std(:,2),'k')
    subplot(613); hold all;
      plot_std(1,gaitCycle,F3_mean_std)
      plot(gaitCycle,F3_mean_std(:,2),'k')
    subplot(614); hold all;
      plot_std(1,gaitCycle,theta_mean_std)
      plot(gaitCycle,theta_mean_std(:,2),'k')
    subplot(615); hold all;
      plot_std(1,gaitCycle,theta_d_mean_std)
      plot(gaitCycle,theta_d_mean_std(:,2),'k')
    subplot(616); hold all;
      plot_std(1,gaitCycle,GRF_mean_std)
      plot(gaitCycle,GRF_mean_std(:,2),'k')


  % --- F1 HS Threshold
  stds = TS.Stats.RMean_std.r_Pa{2}(:,F1_n) ...
          - TS.Stats.RMean_std.r_Pa{3}(:,F1_n);
  F1_HS = max(stds)

  figure;
    subplot(211); hold all;
      title('Mean F1','fontsize',20)
      plot(TS.Stats.RMean_std.r_Pa{2}(:,F1_n))
   subplot(212); hold all;
      title('STD F1','fontsize',20)
      plot(stds)

  % --- F1 HS Threshold
  stds_2 = TS.Stats.RMean_std.r_Pa{2}(:,F2_n) ...
          - TS.Stats.RMean_std.r_Pa{3}(:,F2_n);
  stds_3 = TS.Stats.RMean_std.r_Pa{2}(:,F3_n) ...
          - TS.Stats.RMean_std.r_Pa{3}(:,F3_n);



  figure;
    subplot(211); hold all;
      title('Mean F2 + F3','fontsize',20)
      plot(TS.Stats.RMean_std.r_Pa{2}(:,F2_n) ...
            + TS.Stats.RMean_std.r_Pa{2}(:,F3_n) )
   subplot(212); hold all;
      title('STD F1','fontsize',20)
      plot(stds_2 + stds_3)

    return
end
