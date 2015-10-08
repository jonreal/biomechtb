% Load data files

directory = ['~/Dropbox/Professional/UW_PHD/Prosthetic_Research', ...
            '/Data/PA_A01/Data/'];
dataFiles = [{'PA_SS'}, {'PR_SS'}];
currentVars = whos;

for i=1:numel(dataFiles)
  if ismember(dataFiles{i}, [currentVars(:).name])
    continue;
  end
  load([directory, dataFiles{i}, '/', dataFiles{i}, '.mat'])
end

gaitCycle = linspace(0,100,1001);


palmerDataPath = ['~/Research_link/Prosthetic_Research/Source/', ...
                  'src_v4/data/dataComp/C16Data.mat'];

palmerData = load(palmerDataPath);


t = gaitCycle;
X1_name = 'PA_SS';
X1 = [PA_SS.Stats.RMean_std.RAnkleAngles.X, ...
      PA_SS.Stats.LMean_std.LAnkleAngles.X];

T1 = [PA_SS.Stats.RMean_std.RAnkleMoment.X, ...
      PA_SS.Stats.LMean_std.LAnkleMoment.X];

P1 = [PA_SS.Stats.RMean_std.RAnklePower.X, ...
      PA_SS.Stats.LMean_std.LAnklePower.X];

X2_name = 'PR_SS';
X2 = [PR_SS.Stats.RMean_std.RAnkleAngles.X, ...
      PR_SS.Stats.LMean_std.LAnkleAngles.X];

T2 = [PR_SS.Stats.RMean_std.RAnkleMoment.X, ...
      PR_SS.Stats.LMean_std.LAnkleMoment.X];

P2 = [PR_SS.Stats.RMean_std.RAnklePower.X, ...
      PR_SS.Stats.LMean_std.LAnklePower.X];

X3_name = 'C13';
X3 = [S.R_traj.pos_pos_std, S.R_traj.pos_mean, S.R_traj.pos_neg_std, ...
      S.L_traj.pos_pos_std, S.L_traj.pos_mean, S.L_traj.pos_neg_std];
T3 = -[S.R_traj.mom_pos_std, S.R_traj.mom_mean, S.R_traj.mom_neg_std, ...
      S.L_traj.mom_pos_std, S.L_traj.mom_mean, S.L_traj.mom_neg_std];
P3 = [S.R_traj.pow_pos_std,S.R_traj.pow_mean,S.R_traj.pow_neg_std, ...
      S.L_traj.pow_pos_std,S.L_traj.pow_mean,S.L_traj.pow_neg_std];


h = figure;
xtickFormat = [0 25 50 75 100];

  subplot(331)
    plot_std(h,t,rad2deg(X1));
    title('Powered','interpreter','latex','Fontsize', 22)
    set(gca,'ylim', [-24 48], ...
            'Ytick', [-24 -12 0 12 24], ...
            'Xtick',[0 20 40 60 80 100])
    ylabel({'Angle', '(deg)'},'interpreter','latex','Fontsize',18)
    lh = legend('intact','prosthesis','Location','NorthWest')
    set(lh,'Interpreter','latex', 'Fontsize',15)
    legend boxoff
    grid on
    box on
  subplot(334)
    plot_std(h,t,T1);
    set(gca,'ylim', [-2 1], ...
            'Ytick', [-2 -1 0 1], ...
            'Xtick',[0 20 40 60 80 100])
    ylabel({'Torque','(Nm/kg)'},'interpreter','latex','Fontsize',18)
%    lh = legend('intact','prosthesis','Location','NorthWest')
%    set(lh,'Interpreter','latex', 'Fontsize',13.5)
    legend boxoff
    grid on
    box on
  subplot(337)
    plot_std(h,t,P1);
    set(gca,'ylim', [-2 6], ...
            'Ytick',[-2 0 2 4 6], ...
            'Xtick',[0 20 40 60 80 100])
    ylabel({'Power', '(W/kg)'},'interpreter','latex','Fontsize', 18)
    xlabel('\% Gait','interpreter','latex','Fontsize',18)
%    lh = legend('intact','prosthesis','Location','NorthWest')
%    set(lh,'Interpreter','latex')
    legend boxoff
    grid on
    box on

  subplot(332)
    plot_std(h,t,rad2deg(X2));
    title('Prescribed','interpreter','latex','Fontsize',22)
    set(gca,'ylim', [-24 48], ...
            'Ytick', [-24 -12 0 12 24], ...
            'Xtick',[0 20 40 60 80 100])
    lh = legend('intact','prosthesis','Location','NorthWest')
    set(lh,'Interpreter','latex', 'Fontsize',15)
    legend boxoff
    grid on
    box on
  subplot(335)
    plot_std(h,t,T2);
    set(gca,'ylim',[-2 1], ...
            'Ytick', [-2 -1 0 1], ...
            'Xtick',[0 20 40 60 80 100])
%    lh = legend('intact','prosthesis','Location','NorthWest')
%    set(lh,'Interpreter','latex', 'Fontsize',13.5)
    legend boxoff
    grid on
    box on
  subplot(338)
    plot_std(h,t,P2);
    set(gca,'ylim', [-2 6], ...
            'Ytick',[-2 0 2 4 6], ...
            'Xtick',[0 20 40 60 80 100])
%    lh = legend('intact','prosthesis','Location','NorthWest')
%    set(lh,'Interpreter','latex')
    legend boxoff
    xlabel('\% Gait','interpreter','latex','Fontsize',18)
    grid on
    box on

  subplot(333)
    plot_std(h,t,rad2deg(X3));
    title('Control','interpreter','latex','Fontsize',22)
     set(gca,'ylim', [-24 48], ...
            'Ytick', [-24 -12 0 12 24], ...
            'Xtick',[0 20 40 60 80 100])
    lh = legend('right','left','Location','NorthWest')
    set(lh,'Interpreter','latex','Fontsize',15)
    legend boxoff
    grid on
    box on
  subplot(336)
    plot_std(h,t,T3);
    set(gca,'ylim', [-2 1], ...
            'Ytick', [-2 -1 0 1], ...
            'Xtick',[0 20 40 60 80 100])
%    lh = legend('right','left','Location','NorthWest')
%    set(lh,'Interpreter','latex', 'Fontsize',13.5)
    legend boxoff
    grid on
    box on
  subplot(339)
    plot_std(h,t,P3);
    set(gca,'ylim', [-2 6], ...
            'Ytick',[-2 0 2 4 6], ...
            'Xtick',[0 20 40 60 80 100])
%    lh = legend('right','left','Location','NorthWest')
%    set(lh,'Interpreter','latex')
    legend boxoff
    xlabel('\% Gait','interpreter','latex','Fontsize',18)
    grid on
    box on





