function P = pedar_insole_ui(insole_type_str)
  path = ['~/Dropbox/Professional/UW_PHD/Prosthetic_Research/', ...
          'Source_v2/Biomechanics_Analysis_Toolbox/pedar_insole_ui', ...
          '/insole_type_xy_points'];

  xy_points_file_str = [insole_type_str,'_insole.txt'];
  xy_points = csvread([path,'/',xy_points_file_str]);

  x_points = xy_points(:,1);
  y_points = xy_points(:,2);

  % Allocate Struct
  P.X_points = nan*zeros(200,4);
  P.Y_points = nan*zeros(200,4);
  P.centroid = nan*zeros(200,2);


  X_boundary_points = [];
  Y_boundary_points = [];

  fh = figure; clf;
    plot(x_points,y_points,'ok'), hold all
    title(['Click exterior points.']);

  while(1)
    waitforbuttonpress;
    Cp = get(gca,'CurrentPoint');
    xp = Cp(2,1);
    yp = Cp(2,2);
    [dp,Ip] = min((x_points - xp).^2 + (y_points - yp).^2);
    plot(x_points(Ip),y_points(Ip),'or')

    % Check if point has been clicked
    if ((sum(X_boundary_points == x_points(Ip)) ~= 0) ...
         && (sum(Y_boundary_points == y_points(Ip)) ~= 0))
      break;
    end

    X_boundary_points = [X_boundary_points; x_points(Ip)];
    Y_boundary_points = [Y_boundary_points; y_points(Ip)];
  end

  P.left_boundary = [X_boundary_points, Y_boundary_points];
  P.right_boundary = [-X_boundary_points, Y_boundary_points];

  % 99 cells for each insole
  for i=1:99
    index = nan*ones(1,4);
    while(1)
      figure(fh), clf;
        plot(x_points,y_points,'ok'), hold all
        title(['Click points for cell ', num2str(i)])
      for j=1:i-1
        figure(fh), hold all
          line([P.X_points(j,:),P.X_points(j,1)], ...
               [P.Y_points(j,:),P.Y_points(j,1)], ...
               'Color','k')
          line([P.X_points(j+100,:),P.X_points(j+100,1)], ...
               [P.Y_points(j+100,:),P.Y_points(j+100,1)], ...
               'Color','k')
          plot(P.centroid(j,1),P.centroid(j,2),'or')
          plot(P.centroid(j+100,1),P.centroid(j+100,2),'or')
      end

      for k=1:4
        waitforbuttonpress;

        Cp = get(gca,'CurrentPoint');
        xp = Cp(2,1);
        yp = Cp(2,2);

        [dp,Ip] = min((x_points - xp).^2 + (y_points - yp).^2);
        plot(x_points(Ip),y_points(Ip),'or')
        index(k) = Ip;
      end
      index
      P.X_points(i,:) = x_points(index)
      P.Y_points(i,:) = y_points(index);
      P.centroid(i,:) = [mean(P.X_points(i,:)), mean(P.Y_points(i,:))];

      P.X_points(i+100,:) = -P.X_points(i,:);
      P.Y_points(i+100,:) = P.Y_points(i,:);
      P.centroid(i+100,:) = [mean(-P.X_points(i,:)), mean(P.Y_points(i,:))];

      figure(fh), clf
        plot(x_points,y_points,'ok'), hold all

        line([P.X_points(i,:),P.X_points(i,1)], ...
             [P.Y_points(i,:),P.Y_points(i,1)], ...
             'Color','k')
        line([P.X_points(i+100,:),P.X_points(i+100,1)], ...
             [P.Y_points(i+100,:),P.Y_points(i+100,1)], ...
             'Color','k')
        plot(P.centroid(i,1),P.centroid(i,2),'or')
        plot(P.centroid(i+100,1),P.centroid(i+100,2),'or')
        title('Is this correct')
      varin = input('y/n?', 's');
      if varin == 'y'
        break
      else
        continue
      end
    end
  end
end
