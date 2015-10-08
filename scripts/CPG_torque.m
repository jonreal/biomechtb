% --- Central Pattern Generator - Torque Symmetry
%
%
trials = [{'PA_SS'},{'PA_FS'},{'PA_S'}];
load_trials(trials);

S = PA_SS;

gaitCycle = linspace(0,100,1001);

pros_pos = S.Stats.LMean_std.LAnkleAngles.X;
pros_vel = S.Stats.LMean_std.d_LAnkleAngles.X;
pros_mom = S.Stats.LMean_std.LAnkleMoment.X;

intact_pos = S.Stats.RMean_std.RAnkleAngles.X;
intact_vel = S.Stats.RMean_std.d_RAnkleAngles.X;
intact_mom = S.Stats.RMean_std.RAnkleMoment.X;

TO_pros = S.Stats.LMean_std.TO;
TO_intact = S.Stats.RMean_std.TO;
MDF_pros = S.Stats.LMean_std.MDF;
MDF_intact = S.Stats.RMean_std.MDF;

% --- Piecwise Time Normlized
tt = linspace(0,1,500)';
ii_pros = find(gaitCycle == round(MDF_pros(2)*10)./10);
ii_intact = find(gaitCycle == round(MDF_intact(2)*10)./10);
jj_pros = find(gaitCycle == round(TO_pros(2)*10)./10);
jj_intact = find(gaitCycle == round(TO_intact(2)*10)./10);

gs = gaitCycle(1:ii_pros)';
gc_norm = (gs - gs(1))./(gs(end) - gs(1));
pros_mom_load = pros_mom(1:ii_pros,:);
pros_mom_load = interp1(gc_norm,pros_mom_load,tt,'spline');
pros_pos_load = pros_pos(1:ii_pros,:);
pros_pos_load = interp1(gc_norm,pros_pos_load,tt,'spline');
pros_vel_load = pros_vel(1:ii_pros,:);
pros_vel_load = interp1(gc_norm,pros_vel_load,tt,'spline');

gs = gaitCycle(ii_pros:jj_pros)';
gc_norm = (gs - gs(1))./(gs(end) - gs(1));
pros_mom_unload = pros_mom(ii_pros:jj_pros,:);
pros_mom_unload = interp1(gc_norm,pros_mom_unload,tt,'spline');
pros_pos_unload = pros_pos(ii_pros:jj_pros,:);
pros_pos_unload = interp1(gc_norm,pros_pos_unload,tt,'spline');
pros_vel_unload = pros_vel(ii_pros:jj_pros,:);
pros_vel_unload = interp1(gc_norm,pros_vel_unload,tt,'spline');

gs = gaitCycle(jj_pros:end)';
gc_norm = (gs - gs(1))./(gs(end) - gs(1));
pros_mom_swing = pros_mom(jj_pros:end,:);
pros_mom_swing = interp1(gc_norm,pros_mom_swing,tt,'spline');
pros_pos_swing = pros_pos(jj_pros:end,:);
pros_pos_swing = interp1(gc_norm,pros_pos_swing,tt,'spline');
pros_vel_swing = pros_vel(jj_pros:end,:);
pros_vel_swing = interp1(gc_norm,pros_vel_swing,tt,'spline');

gs = gaitCycle(1:ii_intact)';
gc_norm = (gs - gs(1))./(gs(end) - gs(1));
intact_mom_load = intact_mom(1:ii_intact,:);
intact_mom_load = interp1(gc_norm,intact_mom_load,tt,'spline');
intact_pos_load = intact_pos(1:ii_intact,:);
intact_pos_load = interp1(gc_norm,intact_pos_load,tt,'spline');
intact_vel_load = intact_vel(1:ii_intact,:);
intact_vel_load = interp1(gc_norm,intact_vel_load,tt,'spline');

gs = gaitCycle(ii_intact:jj_intact)';
gc_norm = (gs - gs(1))./(gs(end) - gs(1));
intact_mom_unload = intact_mom(ii_intact:jj_intact,:);
intact_mom_unload = interp1(gc_norm,intact_mom_unload,tt,'spline');
intact_pos_unload = intact_pos(ii_intact:jj_intact,:);
intact_pos_unload = interp1(gc_norm,intact_pos_unload,tt,'spline');
intact_vel_unload = intact_vel(ii_intact:jj_intact,:);
intact_vel_unload = interp1(gc_norm,intact_vel_unload,tt,'spline');

gs = gaitCycle(jj_intact:end)';
gc_norm = (gs - gs(1))./(gs(end) - gs(1));
intact_mom_swing = intact_mom(jj_intact:end,:);
intact_mom_swing = interp1(gc_norm,intact_mom_swing,tt,'spline');
intact_pos_swing = intact_pos(jj_intact:end,:);
intact_pos_swing = interp1(gc_norm,intact_pos_swing,tt,'spline');
intact_vel_swing = intact_vel(jj_intact:end,:);
intact_vel_swing = interp1(gc_norm,intact_vel_swing,tt,'spline');


% -----------------------------------------------------------------------------
% Torque Error (Continous)
% -----------------------------------------------------------------------------
if(1)

  figure; hold all;
    subplot(311); hold all;
      title({'Torque Error','Prosthetic'},'fontsize',20)
      plot_std(gaitCycle,pros_mom,[1,1,1].*0.4)
      plot(gaitCycle,pros_mom(:,2),'k')
      yy = ylim;
      plot([TO_pros(1) TO_pros(1)],[yy(1),yy(2)],'--k')
      plot([TO_pros(2) TO_pros(2)],[yy(1),yy(2)],'-k')
      plot([TO_pros(3) TO_pros(3)],[yy(1),yy(2)],'--k')

      plot([MDF_pros(1) MDF_pros(1)],[yy(1),yy(2)],'--k')
      plot([MDF_pros(2) MDF_pros(2)],[yy(1),yy(2)],'-k')
      plot([MDF_pros(3) MDF_pros(3)],[yy(1),yy(2)],'--k')

      grid on

    subplot(312); hold all;
      title('Intact','fontsize',20)
      plot_std(gaitCycle,intact_mom,[1,1,1].*0.3)
      plot(gaitCycle,intact_mom(:,2),'k')

      yy = ylim;

      plot([TO_intact(1) TO_intact(1)],[yy(1),yy(2)],'--k')
      plot([TO_intact(2) TO_intact(2)],[yy(1),yy(2)],'-k')
      plot([TO_intact(3) TO_intact(3)],[yy(1),yy(2)],'--k')

      plot([MDF_intact(1) MDF_intact(1)],[yy(1),yy(2)],'--k')
      plot([MDF_intact(2) MDF_intact(2)],[yy(1),yy(2)],'-k')
      plot([MDF_intact(3) MDF_intact(3)],[yy(1),yy(2)],'--k')


      grid on
    subplot(313); hold all;
      title('Error Signal','fontsize',20)
      plot_std(gaitCycle,intact_mom - pros_mom,[1 1 1].*.4)
      plot(gaitCycle,intact_mom(:,2) - pros_mom(:,2),'k')
      grid on
end
% -----------------------------------------------------------------------------
% Position Error (Continous)
% -----------------------------------------------------------------------------
if(1)
  figure; hold all;
    subplot(311); hold all;
      title({'Position Error','Prosthetic'},'fontsize',20)
      plot_std(gaitCycle,pros_pos,[1,1,1].*0.4)
      plot(gaitCycle,pros_pos(:,2),'k')
      yy = ylim;
      plot([TO_pros(1) TO_pros(1)],[yy(1),yy(2)],'--k')
      plot([TO_pros(2) TO_pros(2)],[yy(1),yy(2)],'-k')
      plot([TO_pros(3) TO_pros(3)],[yy(1),yy(2)],'--k')

      plot([MDF_pros(1) MDF_pros(1)],[yy(1),yy(2)],'--k')
      plot([MDF_pros(2) MDF_pros(2)],[yy(1),yy(2)],'-k')
      plot([MDF_pros(3) MDF_pros(3)],[yy(1),yy(2)],'--k')

      grid on

    subplot(312); hold all;
      title('Intact','fontsize',20)
      plot_std(gaitCycle,intact_pos,[1,1,1].*0.3)
      plot(gaitCycle,intact_pos(:,2),'k')

      yy = ylim;

      plot([TO_intact(1) TO_intact(1)],[yy(1),yy(2)],'--k')
      plot([TO_intact(2) TO_intact(2)],[yy(1),yy(2)],'-k')
      plot([TO_intact(3) TO_intact(3)],[yy(1),yy(2)],'--k')

      plot([MDF_intact(1) MDF_intact(1)],[yy(1),yy(2)],'--k')
      plot([MDF_intact(2) MDF_intact(2)],[yy(1),yy(2)],'-k')
      plot([MDF_intact(3) MDF_intact(3)],[yy(1),yy(2)],'--k')


      grid on
    subplot(313); hold all;
      title('Error Signal','fontsize',20)
      plot_std(gaitCycle,intact_pos - pros_pos,[1 1 1].*.4)
      plot(gaitCycle,intact_pos(:,2) - pros_pos(:,2),'k')
      grid on
end

% -----------------------------------------------------------------------------
% Velocity Error (Continous)
% -----------------------------------------------------------------------------
if(1)
  figure; hold all;
    subplot(311); hold all;
      title({'Velocity Error','Prosthetic'},'fontsize',20)
      plot_std(gaitCycle,pros_vel,[1,1,1].*0.4)
      plot(gaitCycle,pros_vel(:,2),'k')
      yy = ylim;
      plot([TO_pros(1) TO_pros(1)],[yy(1),yy(2)],'--k')
      plot([TO_pros(2) TO_pros(2)],[yy(1),yy(2)],'-k')
      plot([TO_pros(3) TO_pros(3)],[yy(1),yy(2)],'--k')

      plot([MDF_pros(1) MDF_pros(1)],[yy(1),yy(2)],'--k')
      plot([MDF_pros(2) MDF_pros(2)],[yy(1),yy(2)],'-k')
      plot([MDF_pros(3) MDF_pros(3)],[yy(1),yy(2)],'--k')

      grid on

    subplot(312); hold all;
      title('Intact','fontsize',20)
      plot_std(gaitCycle,intact_vel,[1,1,1].*0.3)
      plot(gaitCycle,intact_vel(:,2),'k')

      yy = ylim;

      plot([TO_intact(1) TO_intact(1)],[yy(1),yy(2)],'--k')
      plot([TO_intact(2) TO_intact(2)],[yy(1),yy(2)],'-k')
      plot([TO_intact(3) TO_intact(3)],[yy(1),yy(2)],'--k')

      plot([MDF_intact(1) MDF_intact(1)],[yy(1),yy(2)],'--k')
      plot([MDF_intact(2) MDF_intact(2)],[yy(1),yy(2)],'-k')
      plot([MDF_intact(3) MDF_intact(3)],[yy(1),yy(2)],'--k')


      grid on
    subplot(313); hold all;
      title('Error Signal','fontsize',20)
      plot_std(gaitCycle,intact_vel - pros_vel,[1 1 1].*.4)
      plot(gaitCycle,intact_vel(:,2) - pros_vel(:,2),'k')
      grid on
end

% -----------------------------------------------------------------------------
% Piece Wise Torque
% -----------------------------------------------------------------------------
if(1)

  figure; hold all
    subplot(231); hold all;
      title('Loading','fontsize',20)
      plot(tt,pros_mom_load(:,2),'-k')
      plot(tt,intact_mom_load(:,2),'-k')
      plot_std(tt',pros_mom_load,[1,0,1].*0.5)
      plot_std(tt',intact_mom_load,[1,1,1].*0.4)
      grid on
      ylim([-2,0.5])

    subplot(232); hold all;
      title({'Piecwise Time Norm Torque','Unloading'},'fontsize',20)
      plot(tt,pros_mom_unload(:,2),'-k')
      plot(tt,intact_mom_unload(:,2),'-k')
      plot_std(tt',pros_mom_unload,[1,0,1].*0.5)
      plot_std(tt',intact_mom_unload,[1,1,1].*0.4)
      grid on

      ylim([-2,0.5])
    subplot(233); hold all;
      title('Swing','fontsize',20)
      plot(tt,pros_mom_swing(:,2),'-k')
      plot(tt,intact_mom_swing(:,2),'--k')
      plot_std(tt',pros_mom_swing,[1,0,1].*0.5)
      plot_std(tt',intact_mom_swing,[1,1,1].*0.4)
      grid on
      ylim([-2,0.5])

    subplot(234); hold all;
      title('Loading','fontsize',20)
      plot(tt,intact_mom_load(:,2)-pros_mom_load(:,2),'k')
      plot_std(tt',intact_mom_load-pros_mom_load,[1,1,1].*0.1)
      grid on
      ylim([-2,0.5])

    subplot(235); hold all;
      title({'Error Signal','Unloading'},'fontsize',20)
      plot(tt,intact_mom_unload(:,2)-pros_mom_unload(:,2),'k')
      plot_std(tt',intact_mom_unload-pros_mom_unload,[1,1,1].*0.1)
      grid on
      ylim([-2,0.5])
    subplot(236); hold all;
      title('Swing','fontsize',20)
      plot(tt,intact_mom_swing(:,2)-pros_mom_swing(:,2),'k')
      plot_std(tt',intact_mom_swing-pros_mom_swing,[1,1,1].*0.1)
      grid on
      ylim([-2,0.5])
end


% -----------------------------------------------------------------------------
% Piece Wise Position
% -----------------------------------------------------------------------------
if(1)

  figure; hold all
    subplot(231); hold all;
      title('Loading','fontsize',20)
      plot(tt,pros_pos_load(:,2),'-k')
      plot(tt,intact_pos_load(:,2),'-k')
      plot_std(tt',pros_pos_load,[1,0,1].*0.5)
      plot_std(tt',intact_pos_load,[1,1,1].*0.4)
      grid on
      ylim([-0.5,0.5])

    subplot(232); hold all;
      title({'Piecwise Time Norm Position','Unloading'},'fontsize',20)
      plot(tt,pros_pos_unload(:,2),'-k')
      plot(tt,intact_pos_unload(:,2),'-k')
      plot_std(tt',pros_pos_unload,[1,0,1].*0.5)
      plot_std(tt',intact_pos_unload,[1,1,1].*0.4)
      grid on
      ylim([-0.5,0.5])

    subplot(233); hold all;
      title('Swing','fontsize',20)
      plot(tt,pros_pos_swing(:,2),'-k')
      plot(tt,intact_pos_swing(:,2),'--k')
      plot_std(tt',pros_pos_swing,[1,0,1].*0.5)
      plot_std(tt',intact_pos_swing,[1,1,1].*0.4)
      grid on
      ylim([-0.5,0.5])

    subplot(234); hold all;
      title('Loading','fontsize',20)
      plot(tt,intact_pos_load(:,2)-pros_pos_load(:,2),'k')
      plot_std(tt',intact_pos_load-pros_pos_load,[1,1,1].*0.1)
      grid on
      ylim([-0.5,0.5])

    subplot(235); hold all;
      title({'Error Signal','Unloading'},'fontsize',20)
      plot(tt,intact_pos_unload(:,2)-pros_pos_unload(:,2),'k')
      plot_std(tt',intact_pos_unload-pros_pos_unload,[1,1,1].*0.1)
      grid on
      ylim([-0.5,0.5])
    subplot(236); hold all;
      title('Swing','fontsize',20)
      plot(tt,intact_pos_swing(:,2)-pros_pos_swing(:,2),'k')
      plot_std(tt',intact_pos_swing-pros_pos_swing,[1,1,1].*0.1)
      grid on
      ylim([-0.5,0.5])
end

% -----------------------------------------------------------------------------
% Piece Wise Velocity
% -----------------------------------------------------------------------------
if(1)

  figure; hold all
    subplot(231); hold all;
      title('Loading','fontsize',20)
      plot(tt,pros_vel_load(:,2),'-k')
      plot(tt,intact_vel_load(:,2),'-k')
      plot_std(tt',pros_vel_load,[1,0,1].*0.5)
      plot_std(tt',intact_vel_load,[1,1,1].*0.4)
      grid on
      ylim([-6,4])

    subplot(232); hold all;
      title({'Piecwise Time Norm Velocity','Unloading'},'fontsize',20)
      plot(tt,pros_vel_unload(:,2),'-k')
      plot(tt,intact_vel_unload(:,2),'-k')
      plot_std(tt',pros_vel_unload,[1,0,1].*0.5)
      plot_std(tt',intact_vel_unload,[1,1,1].*0.4)
      grid on
      ylim([-6,4])

    subplot(233); hold all;
      title('Swing','fontsize',20)
      plot(tt,pros_vel_swing(:,2),'-k')
      plot(tt,intact_vel_swing(:,2),'--k')
      plot_std(tt',pros_vel_swing,[1,0,1].*0.5)
      plot_std(tt',intact_vel_swing,[1,1,1].*0.4)
      grid on
      ylim([-6,4])

    subplot(234); hold all;
      title('Loading','fontsize',20)
      plot(tt,intact_vel_load(:,2)-pros_vel_load(:,2),'k')
      plot_std(tt',intact_vel_load-pros_vel_load,[1,1,1].*0.1)
      grid on
      ylim([-6,4])

    subplot(235); hold all;
      title({'Error Signal','Unloading'},'fontsize',20)
      plot(tt,intact_vel_unload(:,2)-pros_vel_unload(:,2),'k')
      plot_std(tt',intact_vel_unload-pros_vel_unload,[1,1,1].*0.1)
      grid on
      ylim([-6,4])
    subplot(236); hold all;
      title('Swing','fontsize',20)
      plot(tt,intact_vel_swing(:,2)-pros_vel_swing(:,2),'k')
      plot_std(tt',intact_vel_swing-pros_vel_swing,[1,1,1].*0.1)
      grid on
      ylim([-6,4])
end

% -----------------------------------------------------------------------------
% Stiffness
% -----------------------------------------------------------------------------
if(0)
  figure; hold all;
    title('Stiffness','fontsize',20)
    plot(pros_pos(:,2),pros_mom(:,2),'-r')
    plot(intact_pos(:,2),intact_mom(:,2),'-k')
    grid on
    xlim([-0.3,0.4]);
    ylim([-2,0.5]);

  figure; hold all;
    title('Torque-Velocity','fontsize',20)
    plot(pros_vel(:,2),pros_mom(:,2),'-r')
    plot(intact_vel(:,2),intact_mom(:,2),'-k')
    grid on

  figure;
    subplot(131); hold all;
      plot(pros_pos_load(:,2),pros_mom_load(:,2),'-r')
      plot(intact_pos_load(:,2),intact_mom_load(:,2),'-k')
      grid on
      xlim([-0.3,0.4]);
      ylim([-2,0.5]);

    subplot(132); hold all;
      title('Stiffness','fontsize',20)
      plot(pros_pos_unload(:,2),pros_mom_unload(:,2),'-r')
      plot(intact_pos_unload(:,2),intact_mom_unload(:,2),'-k')
      grid on
      xlim([-0.3,0.4]);
      ylim([-2,0.5]);
    subplot(133); hold all;
      plot(pros_pos_swing(:,2),pros_mom_swing(:,2),'-r')
      plot(intact_pos_swing(:,2),intact_mom_swing(:,2),'-k')
      grid on
      xlim([-0.3,0.4]);
      ylim([-2,0.5]);

  figure;
    subplot(131); hold all;
      plot(pros_vel_load(:,2),pros_mom_load(:,2),'-r')
      plot(intact_vel_load(:,2),intact_mom_load(:,2),'-k')
      grid on
    subplot(132); hold all;
      title('Torque-Velocity','fontsize',20)
      plot(pros_vel_unload(:,2),pros_mom_unload(:,2),'-r')
      plot(intact_vel_unload(:,2),intact_mom_unload(:,2),'-k')
      grid on
    subplot(133); hold all;
      plot(pros_vel_swing(:,2),pros_mom_swing(:,2),'-r')
      plot(intact_vel_swing(:,2),intact_mom_swing(:,2),'-k')
      grid on
end

