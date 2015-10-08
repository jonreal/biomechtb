function plot_trial_ankle_knee_hip(trial)
  directory = ['~/Dropbox/Professional/UW_PHD/Prosthetic_Research', ...
               '/Data/PA_A01/Data/'];

  S = load([directory, trial, '/', trial, '.mat']);
  t = linspace(0,100,1001);

%  S.([genvarname(trial)]).Stats.RMean_std

  % Allocate Data
  X1 = [S.([genvarname(trial)]).Stats.RMean_std.RAnkleAngles.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LAnkleAngles.X];
  X2 = [S.([genvarname(trial)]).Stats.RMean_std.RKneeAngles.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LKneeAngles.X];
  X3 = [S.([genvarname(trial)]).Stats.RMean_std.RHipAngles.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LHipAngles.X];

  T1 = [S.([genvarname(trial)]).Stats.RMean_std.RAnkleMoment.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LAnkleMoment.X];
  T2 = [S.([genvarname(trial)]).Stats.RMean_std.RKneeMoment.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LKneeMoment.X];
  T3 = [S.([genvarname(trial)]).Stats.RMean_std.RHipMoment.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LHipMoment.X];

  P1 = [S.([genvarname(trial)]).Stats.RMean_std.RAnklePower.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LAnklePower.X];
  P2 = [S.([genvarname(trial)]).Stats.RMean_std.RKneePower.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LKneePower.X];
  P3 = [S.([genvarname(trial)]).Stats.RMean_std.RHipPower.X, ...
        S.([genvarname(trial)]).Stats.LMean_std.LHipPower.X];


  % Plot
  h = figure;
    subplot(331)
    plot_std(h,t,rad2deg(X1));
    title('Ankle','interpreter','latex','fontsize',20)
    set(gca,'ylim', [-24 24], ...
            'Ytick', [-24 -12 0 12 24], ...
            'Xtick',[0 20 40 60 80 100])
    ylabel('Angle (deg)','interpreter','latex','Fontsize',14)
    lh = legend('intact','prosthesis','Location','SouthWest');
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

  subplot(334)
    plot_std(h,t,T1);
    set(gca,'ylim', [-2 1], ...
            'Ytick', [-2 -1 0 1], ...
            'Xtick',[0 20 40 60 80 100])
    ylabel('Torque (Nm/kg)','interpreter','latex','Fontsize',14)
    lh = legend('intact','prosthesis','Location','SouthWest');
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

  subplot(337)
    plot_std(h,t,P1);
    set(gca,'ylim', [-2 6], ...
            'Ytick',[-2 0 2 4 6], ...
            'Xtick',[0 20 40 60 80 100])
    ylabel('Power (W/kg)','interpreter','latex','Fontsize', 14)
    xlabel('\% Gait','interpreter','latex','Fontsize',14)
    lh = legend('intact','prosthesis','Location','NorthWest');
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on


  subplot(332)
    plot_std(h,t,rad2deg(X2));
    title('Knee','interpreter','latex','Fontsize',20)
    set(gca,'ylim', [-10 80], ...
            'Ytick', [-10 -0 20 40 60 80], ...
            'Xtick',[0 20 40 60 80 100])
    lh = legend('intact','prosthesis','Location','NorthWest')
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

  subplot(335)
    plot_std(h,t,T2);
    set(gca,'ylim', [-2 1], ...
            'Ytick', [-2 -1 0 1], ...
            'Xtick',[0 20 40 60 80 100])
    lh = legend('intact','prosthesis','Location','SouthWest')
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

  subplot(338)
    plot_std(h,t,P2);
    set(gca,'ylim', [-2 6], ...
            'Ytick',[-2 0 2 4 6], ...
            'Xtick',[0 20 40 60 80 100])
    xlabel('\% Gait','interpreter','latex','Fontsize',14)
    lh = legend('intact','prosthesis','Location','NorthWest')
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on


  subplot(333)
    plot_std(h,t,rad2deg(X3));
    title('Hip','interpreter','latex','Fontsize',20)
     set(gca,'ylim', [-20 60], ...
            'Ytick', [-20 0 20 40 60], ...
            'Xtick',[0 20 40 60 80 100])
    lh = legend('intact','prosthesis','Location','North')
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

  subplot(336)
    plot_std(h,t,T3);
    set(gca,'ylim', [-2 1], ...
            'Ytick', [-2 -1 0 1], ...
            'Xtick',[0 20 40 60 80 100])
    lh = legend('intact','prosthesis','Location','SouthWest')
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

  subplot(339)
    plot_std(h,t,P3);
    set(gca,'ylim', [-2 6], ...
            'Ytick',[-2 0 2 4 6], ...
            'Xtick',[0 20 40 60 80 100])
    xlabel('\% Gait','interpreter','latex','Fontsize',14)
    lh = legend('intact','prosthesis','Location','NorthWest')
    set(lh,'Interpreter','latex')
    legend boxoff, grid on, box on

    suph = suptitle('PR-SL')
    set(suph, 'interpreter','latex','Fontsize',20)
end
