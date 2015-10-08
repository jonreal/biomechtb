% --- Neural Net Scripts
%
%

dataFile = [{'PA_SS'},{'PA_FS'},{'PA_SL_01'}, ...
             {'PA_C'},{'PA_S'},{'PA_TSS_2'}, ...
            {'PR_SS'},{'PR_SL'}];

load_trials(dataFile);

% -----------------------------------------------------------------------------
% Inputs
% -----------------------------------------------------------------------------

% Sensor Location
F1_n = 3;
F2_n = 78;
F3_n = 98;

% Training, and Validation
TS = {PR_SS};
VS = {PR_SL};

plotSensLoc = 1;

trainNN = 1;
validateNN = 1;
createNetwork = 1;
plotTrain = 1;
plotValidate = 1;

tapDelay = 1:10;
numHiddenLay = 3;

% -----------------------------------------------------------------------------
% Create Network
% -----------------------------------------------------------------------------
if createNetwork
  if numel(tapDelay) > 1
    lrn_net = feedforwardnet(numHiddenLay);
    B
  else
    lrn_net = layrecnet(tapDelay,numHiddenLay);
  end


  lrn_net.layers{1}.transferFcn = 'logsig';
  lrn_net.layers{2}.transferFcn = 'purelin';

  lrn_net.performFcn = 'mse';
  lrn_net.performParam.regularization = 0.0;

  lrn_net.trainFcn = 'trainlm';
  lrn_net.trainParam.epochs = 10000;
  lrn_net.trainParam.min_grad = 0;
  lrn_net.trainParam.mu_max = 1e100;
  lrn_net.trainParam.showCommandLine = true;

  lrn_net.inputs{1}.processFcns = {};
  lrn_net.outputs{2}.processFcns = {};
end

if(trainNN)
  % ---------------------------------------------------------------------------
  % Training/Validation Set
  % ---------------------------------------------------------------------------
  D = {TS,VS};
  X = {{},{}};
  y = {{},{}};
  for i=1:2
    for k=1:numel(D{i})
      r_gaitPhase = D{i}{k}.Events.R.gaitPhase;
      r_gaitPhase = r_gaitPhase.*(r_gaitPhase == 1);
      X{i}{k} = D{i}{k}.Pedar.Interp.r_Pa(:,[F1_n,F2_n,F3_n]);
     % y{i}{k} = [r_gaitPhase];
      y{i}{k} = D{i}{k}.Model.RAnkleMoment.X;
    end
  end

  % ---------------------------------------------------------------------------
  % Training
  % ---------------------------------------------------------------------------
  u_cat = con2seq(X{1}{1}');
  t_cat = con2seq(y{1}{1}');

  for k=2:numel(X{1});
    u_cat = catsamples(u_cat,con2seq(X{1}{k}'),'pad');
    t_cat = catsamples(t_cat,con2seq(y{1}{k}'),'pad');
  end

  [p,Pi,Ai,t] = preparets(lrn_net,u_cat,t_cat);
  lrn_net = train(lrn_net,p,t,Pi);

  yp = cell(1,numel(X{1}));
  for k=1:numel(X{1})
    o = lrn_net(con2seq(X{1}{k}'));
    yp{k} = cell2mat(o)';
  end
end


% -----------------------------------------------------------------------------
% Testing DTRNN
% -----------------------------------------------------------------------------
if(validateNN)
  yv = cell(1,numel(X{2})); 
  for k=1:numel(X{2});
    u = con2seq(X{2}{k}');
    o = lrn_net(u);
    yv{k} = cell2mat(o)';
  end
end

% -----------------------------------------------------------------------------
% Plots
% -----------------------------------------------------------------------------
if (plotTrain)
  for k=1:numel(D{1})
    [~,n] = size(y{1}{1});
    str = strsplit(D{1}{k}.name,'_');
    figure;
      subplot(2*n,1,1:n),
      plot(X{1}{k})
      title(['Traing Set: ',str{1},'-',str{2},'- Inputs'],'fontsize',20)
      xlim([0,numel(y{1}{k})]);
      axis off;

      subplot(2*n,1,n+1); hold all
      title('Ouputs','fontsize',20)
      plot(y{1}{k}(:,1),'k')
      plot(yp{k}(:,1),'--r'),
      xlim([0,numel(y{1}{k}(:,1))]);
      axis off
      grid on
    for i=2:n
      subplot(2*n,1,n+i); hold all
        plot(y{1}{k}(:,i),'k')
        plot(yp{k}(:,i),'--r'),
        xlim([0,numel(y{1}{k}(:,1))]);
        axis off
    end
  end
end

if (plotValidate)
  for k=1:numel(D{2})
    [~,n] = size(y{2}{1});
    str = strsplit(D{2}{k}.name,'_');
    figure;
      subplot(2*n,1,1:n);
      plot(X{2}{k})
      title(['Validation Set: ',str{1},'-',str{2},' - Inputs'],'fontsize',20)
      xlim([0,numel(y{2}{k})]);
      axis off;

      subplot(2*n,1,n+1); hold all
      title('Ouputs','fontsize',20)
      plot(y{2}{k}(:,1),'k')
      plot(yv{k}(:,1),'--r'),
      xlim([0,numel(y{2}{k}(:,1))]);
      axis off
      grid on
    for i=2:n
      subplot(2*n,1,n+i); hold all
        plot(y{2}{k}(:,i),'k'),
        plot(yv{k}(:,i),'--r'),
        xlim([0,numel(y{2}{k}(:,1))]);
        axis off
    end
  end
end


% -----------------------------------------------------------------------------
% Plot sensor locations
% -----------------------------------------------------------------------------

if(plotSensLoc)
tt = linspace(0,2*pi,100);
circle_x =@(x0) 3.5*cos(tt) + x0;
circle_y =@(y0) 3.5*sin(tt) + y0;

figure, hold all;
  line([PA_SS.Pedar.sensorTemplate.left_boundary(:,1); ...
        PA_SS.Pedar.sensorTemplate.left_boundary(1,1)], ...
       [PA_SS.Pedar.sensorTemplate.left_boundary(:,2); ...
        PA_SS.Pedar.sensorTemplate.left_boundary(1,2)],'Color','k');
  line([PA_SS.Pedar.sensorTemplate.right_boundary(:,1); ...
        PA_SS.Pedar.sensorTemplate.right_boundary(1,1)], ...
       [PA_SS.Pedar.sensorTemplate.right_boundary(:,2); ...
        PA_SS.Pedar.sensorTemplate.right_boundary(1,2)],'Color','k');

  for i=1:numel(PA_SS.Pedar.Interp.l_Pa(1,:))
    patch(circle_x(PA_SS.Pedar.sensorTemplate.l_x_centroid(i)), ...
                   circle_y(PA_SS.Pedar.sensorTemplate.l_y_centroid(i)), ...
                   0);
    patch(circle_x(PA_SS.Pedar.sensorTemplate.r_x_centroid(i)), ...
                   circle_y(PA_SS.Pedar.sensorTemplate.r_y_centroid(i)), ...
                   0);
  end
  patch(circle_x(PA_SS.Pedar.sensorTemplate.r_x_centroid(F1_n)), ...
                 circle_y(PA_SS.Pedar.sensorTemplate.r_y_centroid(F1_n)), ...
                 1);
  patch(circle_x(PA_SS.Pedar.sensorTemplate.r_x_centroid(F2_n)), ...
                 circle_y(PA_SS.Pedar.sensorTemplate.r_y_centroid(F2_n)), ...
                 1);
  patch(circle_x(PA_SS.Pedar.sensorTemplate.r_x_centroid(F3_n)), ...
                 circle_y(PA_SS.Pedar.sensorTemplate.r_y_centroid(F3_n)), ...
                 1);

  axis equal
  axis off
  box off

  return

end

