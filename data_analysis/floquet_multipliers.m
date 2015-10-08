function [maxFM,FM] = floquet_multipliers(t,Sn,numOfMaps,periodicEventLocs)
  debug = 0;

  if debug
    f1 = figure;
    f2 = figure;
  end

  % k-datapoints, m-states
  [k,m] = size(Sn);

  numOfPeriods = numel(periodicEventLocs)-1;
  gaitCycle = linspace(0,100,numOfMaps);
  Sn_norm = zeros(numOfMaps,m,numOfPeriods);
  FM = zeros(numOfMaps,m);
  maxFM = zeros(1,numOfMaps);

  % Interpolate each stride
  for i=1:numOfPeriods
    sIndx = periodicEventLocs(i);
    eIndx = periodicEventLocs(i+1);

    rawTime = t(sIndx:eIndx);
    rawGaitCycle = (rawTime - rawTime(1))./(rawTime(end) - rawTime(1))*100;

    Sn_norm(:,:,i) = ...
      interp1(rawGaitCycle,Sn(sIndx:eIndx,:),gaitCycle,'spline');

    if debug
      figure(f1); clf; hold all;
        plot(gaitCycle,Sn_norm(:,1,i))
        xlabel('% gait')
        ylabel('state')
        title('Time Normalized Trajectories')
        pause(0.1)

      figure(f2); clf; hold all;
        plot3(Sn_norm(:,1,i),Sn_norm(:,2,i),Sn_norm(:,3,i))
        xlabel('Sn_1')
        ylabel('Sn_2')
        zlabel('Sn_3')
        title('Time Normalized Trajectories')
        view(3)
        pause(0.1)

    end
  end

  % Find Fixed Points
  Sn_star = mean(Sn_norm,3);

  if debug
   figure(f1); clf; hold all;
      plot(gaitCycle,squeeze(Sn_norm(:,1,:)))
      plot(gaitCycle,Sn_star(:,1),'--k','LineWidth',3)
      xlabel('% gait')
      ylabel('state')
      title('Fixed Points')
      pause(0.1)

   figure(f2); clf; hold all;
      plot3(squeeze(Sn_norm(:,1,:)),squeeze(Sn_norm(:,2,:)), ...
           squeeze(Sn_norm(:,3,:)))
      plot3(Sn_star(:,1),Sn_star(:,2),Sn_star(:,3),'--k','LineWidth',3)
      xlabel('Sn_1')
      ylabel('Sn_2')
      zlabel('Sn_3')
      title('Fixed Points')
      view(3)
      pause(0.1)

  end

  % Least Squares each map
  for i=1:numOfMaps
    Sk = squeeze(Sn_norm(i,:,:))';

    y = Sk(2:end,:) - Sn_star(i*ones(numOfPeriods-1,1),:);
    y = y(:);

    A = kron(eye(m),Sk(1:end-1,:) - Sn_star(i*ones(numOfPeriods-1,1),:));

    % Solve
    if (rank(A'*A) < numOfPeriods)
      x_opt = pinv(A'*A)*A'*y;
    else
      x_opt = (A'*A)\A'*y;
    end
    FM(i,:) = eig(reshape(x_opt,m,m));
    maxFM(i) = max(abs(FM(i,:)));

    if debug
      figure(f1); clf; hold all;
        plot(gaitCycle,squeeze(Sn_norm(:,1,:)))
        plot(gaitCycle,Sn_star(:,1),'--k','LineWidth',3)
        plot(gaitCycle(i),Sk(:,1)','or')
        pause(0.001)

      figure(f2); clf; hold all;
        plot3(squeeze(Sn_norm(:,1,:)),squeeze(Sn_norm(:,2,:)), ...
              squeeze(Sn_norm(:,3,:)))
        plot3(Sn_star(:,1),Sn_star(:,2),Sn_star(:,3),'--k','LineWidth',3)
        plot3(Sk(:,1)',Sk(:,2)',Sk(:,3)','or')
        view(3);
        pause(0.001)

    end
end
