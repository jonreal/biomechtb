function plot_cop(trialName, sensor_type)
  sensorTemplateFilePath = ['~/Research_link/Prosthetic_Research/', ...
              'Source_v2/Biomechanics_Analysis_Toolbox/', ...
              'pedar_insole_ui/insole_mat_files/', ...
              sensor_type,'_insole.mat'];

  load(sensorTemplateFilePath);

  P

  l_x_centroid = P.centroid(1:99,1)' + 75;
  l_y_centroid = P.centroid(1:99,2)';

  r_x_centroid = P.centroid((1:99) + 100,1)' - 75;
  r_y_centroid = P.centroid((1:99) + 100,2)';

  xycop_func =@(x,y,Pa) [sum(x.*Pa)/sum(Pa), sum(y.*Pa)/sum(Pa)];

  dataFilePath = ['~/Research_link/Prosthetic_Research/',...
                  'Data/PA_A01/Pedar_RAW/', trialName,'.asc'];

  RAW = dlmread(dataFilePath,'\t',7,0);
  time = RAW(:,1);
  l_Pa = RAW(:,2:100);
  r_Pa = RAW(:,101:end-1);
  l_xy_cop = nan*zeros(length(l_Pa),2);
  r_xy_cop = nan*zeros(length(r_Pa),2);

  l_Pa = l_Pa./max(max(l_Pa));
  r_Pa = r_Pa./max(max(r_Pa));

  l_thres = 0.2;
  r_thres = 0.2;
  l_Pa = l_Pa - l_thres;
  l_Pa(l_Pa < 0) = 0;
  r_Pa = r_Pa - r_thres;
  r_Pa(r_Pa < 0) = 0;

  if(0)
  % Plots for threshold
  figure, hold all
    title('Left')
    plot(time,l_Pa)
    plot(time,l_thres*ones(size(time)),'k');
  figure, hold all
    title('Left')
    plot(time,l_Pa)

  figure, hold all
    title('right')
    plot(time,r_Pa)
    plot(time,r_thres*ones(size(time)),'k');
  figure, hold all
    title('right')
    plot(time,r_Pa)
    return
  end

  for i=1:numel(l_Pa(:,1))
  % Calcluate COP from centroid data
    l_xy_cop(i,:) = xycop_func(l_x_centroid,l_y_centroid,l_Pa(i,:));
    r_xy_cop(i,:) = xycop_func(r_x_centroid,r_y_centroid,r_Pa(i,:));
  end

  % Position
figure, hold all
  title('Position','Fontsize', 20)
  plot(time,-l_xy_cop(:,1),'k')
  plot(time,r_xy_cop(:,1),'r')
  xlabel('Time (s)','Fontsize', 20)
  ylabel('COP_x (mm)','Fontsize', 20)
  lh = legend('Left Foot', 'Right Foot')
  set(lh,'FontSize',14);
  set(gca, 'FontSize',14)
  box on
  grid on

figure, hold all
  title('Position','Fontsize', 20)
  plot(time,l_xy_cop(:,2),'k')
  plot(time,r_xy_cop(:,2),'r')
  xlabel('Time (s)','Fontsize', 20)
  ylabel('COP_y (mm)','Fontsize', 20)
  lh = legend('Left Foot', 'Right Foot')
  set(lh,'FontSize',14);
  set(gca, 'FontSize',14)
  grid on
  box on
  ylim([0, 250])

figure, hold all
  % Plot COP
  title('Position', 'FontSize', 20)
  plot(l_xy_cop(:,1), l_xy_cop(:,2), 'k')
  plot(r_xy_cop(:,1), r_xy_cop(:,2), 'r')

%  plot(-l_xy_cop(:,1), l_xy_cop(:,2), 'r--')
%  plot(-r_xy_cop(:,1), r_xy_cop(:,2), 'g--')


  P.left_boundary = [P.left_boundary(:,1) + 75, P.left_boundary(:,2)]
  P.right_boundary = [P.right_boundary(:,1) - 75, P.right_boundary(:,2)]

  % Plot polygon outlines
  left_poly_plot = [P.left_boundary; P.left_boundary(1,:)];
  right_poly_plot = [P.right_boundary; P.right_boundary(1,:)];
  line(left_poly_plot(:,1),left_poly_plot(:,2),'Color','k')
  line(right_poly_plot(:,1),right_poly_plot(:,2),'Color','k')
  xlabel('x (mm)','Fontsize', 20)
  ylabel('y (mm)','Fontsize', 20)
  axis equal
  set(gca, 'FontSize',14)
  grid on

%  set(gca, 'XTick', []);
%  set(gca, 'YTick', []);


end
