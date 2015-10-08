matFilePath = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/', ...
               'AnkleOS/BB_stiffness_switching.mat'];

load(matFilePath);


k_s = 172000;
k_e = 186000;

k_s = 164350;
k_e = 178000;

time = S.time(k_s:k_e) - S.time(k_s);
ankle_pos = S.ankle_pos(k_s:k_e);
cmd = S.Kp(k_s:k_e);
offset = 14;
numMaps = 100;

[pks, locs] = findpeaks(ankle_pos);

figure, hold all
  subplot(211)
  plot(time,cmd, 'k')
  xlabel('Sample','interpreter', ...
         'latex','fontsize',15)
  ylabel('Propotional Gain, $K_p$','interpreter','latex', ...
         'fontsize',15)
  box on
  grid on

  indx = find(pks>-0.02);
  locs = locs(indx);
  pks = pks(indx);

  for i=2:numel(locs)
    if (mod(i,2) == 0)
      ankle_pos(locs(i-1):locs(i)) = -ankle_pos(locs(i-1):locs(i));
    end
  end


  subplot(212), hold on
  plot(time, rad2deg(ankle_pos) + offset, 'k')


  [pks, locs] = findpeaks(rad2deg(ankle_pos) + offset);
  indx = find(pks>15);
  locs = locs(indx);
  pks = pks(indx);
  plot(time(locs),rad2deg(ankle_pos(locs)) + offset,'or')


  T = time(locs(2:end)) - time(locs(1:end-1));
  t0 = time(locs(1));

  % Find peaks that are close together
  indx = find(locs(2:end) - locs(1:end-1) < 200)
  locs(indx) = [];

  xlabel('Time (s)','interpreter', ...
         'latex','fontsize',15)
  ylabel('Angle Angle (deg)','interpreter','latex', ...
         'fontsize',15)
  box on
  grid on

return

if(0)
%FFT
w = my_fft(x,1,1,1);
meanT = ceil(1/w);

tau = 50;
m = 3;

Sn = delay_reconstruct(x,tau,m);
d = nndist(Sn,meanT)

return
end

%% Period by Period - fouier fit
t = time;
x  = ankle_pos + deg2rad(offset);


fh = figure, hold all
  plot(t,x,'k')
  plot(t(locs),x(locs),'or')

T_est = zeros(numel(locs)-1,1);
coeffs = zeros(18,numel(locs)-1);

for i=1:(numel(locs)-1)
  t_win = t(locs(i):locs(i+1));
  x_win = x(locs(i):locs(i+1));
  [f,gof] = fit(t_win,x_win,'fourier8');


  t_est(i) = t(locs(i+1));
  T_est(i) = abs(t(locs(i)) - t(locs(i+1)));
  coeffs(:,i) = coeffvalues(f)';

  figure(fh), clf, hold all
    plot(t,x,'k')
    plot(t(locs(i):locs(i+1)), ...
         x(locs(i):locs(i+1)),'--r')
    plot(t(locs(i:i+1)), ...
         x(locs(i:i+1)),'or')

  pause(0.001)

  if(0)
    figure(fh1), clf,
      plot(f,t_win,x_win)
      pause(0.01)
  end
end

figure;
  plot(T_est);
figure;
  plotyy(t_est,coeffs',t,cmd)
  xlabel('Time (s)','interpreter','latex','fontsize',15)
  ylabel('Coeff','interpreter','latex','fontsize',15)
  title('Windowed Fourier Series', ...
        'interpreter','latex','fontsize',25)


%%
fundFreq = my_fft(ankle_pos,0.005,1,1)



% RLS
t = time;
x  = ankle_pos + deg2rad(offset);


phi_fun =@(t,w) [1, cos(t*w), sin(t*w), ...
                cos(2*t*w), sin(2*t*w), ...
                cos(3*t*w), sin(3*t*w), ...
                cos(4*t*w), sin(4*t*w), ...
                cos(5*t*w), sin(5*t*w), ...
                cos(6*t*w), sin(6*t*w), ...
                cos(7*t*w), sin(7*t*w), ...
                cos(8*t*w), sin(8*t*w)]';


theta0 = zeros(17,1);
P0 = 1000*eye(17,17);
lam = 0.98;

theta_vector = zeros(numel(t),17);
yHat_vector = zeros(numel(t),1);
P_norm_vector = zeros(numel(t),1);
fundFreq_vector = zeros(numel(t),1);

window = 325;
theta_window = zeros(17,5);
covariance_norm_window = zeros(10,1);
sample_window = zeros(window,1);
time_window = zeros(window,1);


fh1 = figure;
t_coeffs = [];
coeffs = [];
gof_coeffs = [];

coeffs(:,end+1) = zeros(1,18);
t_coeffs(end+1) = 0;
gof_coeffs(end+1) = 0;

cnt = 0;
for i=1:numel(t)

  i/numel(t)

  sample_window = [sample_window(2:end); x(i)];
  time_window = [time_window(2:end);, t(i)];

  if ((i>window) && (mod(i,1000) == 0))

    if(cnt == 0)

    [f,gof] = fit(time_window,sample_window, ...
                  'fourier8');
    else

    [f,gof] = fit(time_window,sample_window, ...
            'fourier8','StartPoint', coeffs(:,end));
    end
    cnt = cnt + 1;

    t_coeffs(end+1) = t(i);
    coeffs(:,end+1) = coeffvalues(f)';
    gof_coeffs(end+1) = gof.sse;


    if(1)
    figure(fh1), clf,
%      plot(t(1:i),coeffs(:,1:i)')
      plot(f,time_window,sample_window)
      pause(0.001)
    end

    if(0)
    [pks,loc] = findpeaks(sample_window);

    if(0)
    [pks, loc]
    figure(fh1), clf, hold all
      plot(time_window,rad2deg(sample_window),'k')
      plot(time_window(loc),rad2deg(sample_window(loc)),'or')
      pause(0.001)
    end

    % Threshold
    indx = find(rad2deg(pks)>18 & rad2deg(pks)<22);
    pks = pks(indx);
    loc = loc(indx);

    if(0)
    [pks, loc]
    figure(fh1), clf, hold all
      plot(time_window,rad2deg(sample_window),'k')
      plot(time_window(loc),rad2deg(sample_window(loc)),'or')
      pause(0.001)
    end

    % Elimate double peaks
    indx = find(abs(loc(2:end)-loc(1:end-1)) < 100);
    pks(indx) = [];
    loc(indx) = [];

    if(0)
    [pks, loc]
    figure(fh1), clf, hold all
      plot(time_window,rad2deg(sample_window),'k')
      plot(time_window(loc),rad2deg(sample_window(loc)),'or')
      pause(0.001)
    end

    % Sort in Time
    [loc, indx] = sort(loc);
    pks = pks(indx);

    if(0)
    [pks, loc]
    figure(fh1), clf, hold all
      plot(time_window,rad2deg(sample_window),'k')
      plot(time_window(loc(end-1:end)), ...
           rad2deg(sample_window(loc(end-1:end))),'*g')
      pause(0.25)
    end

    T_est = abs(time_window(loc(end)) - time_window(loc(end-1)));
    fundFreq = 1/T_est;
    fundFreq_vector(i) = fundFreq;

    if (T_est < 1)
      T_est = T_est
      disp('T_est < 1')
      pause
    end
  end
  else
    fundFreq_vector(i) = fundFreq;
  end


  if(1)
  phi = phi_fun(t(i),fundFreq);
  [yHat, thetaHat, P] = rls(17,theta0,P0,lam,x(i),phi);

  theta0 = thetaHat;
  P0 = P;

  P_norm_vector(i) = norm(P,2);

  theta_vector(i,:) = theta0;
  yHat_vector(i) = yHat;

  theta_window = [theta_window(:,2:end),theta0];
  covariance_norm_window = [covariance_norm_window(2:end); norm(P,2)];

  covar_vector(i) = var(covariance_norm_window);
  end
  if(0)
    figure(fh1),
      clf, hold all
      plot(t(1:i),x(1:i), 'k')
      plot(time_window,sample_window,'--r')

%  if( i > 2000)
%
%  figure(fh1);
    %semilogy(var(theta_window')')
    %ylim([1e-10,1e5])
%    plot(t(i-2000:i),P_norm_vector(i-2000:i))
%    title(sprintf('%0.2f',t(i)))
%    pause(0.001)
  end
end

figure;
  plotyy(t_coeffs,coeffs,t,cmd)

figure;
  plotyy(t,P_norm_vector,t,cmd)

figure;
  plotyy(t,theta_vector',t,cmd)

figure;
    plot(t,x,'-k'), hold on
    plot(t,yHat_vector,'--r')

figure;
  plotyy(t,(x - yHat_vector).^2,t,cmd)


return

figure;
  plot(x(:,1),x(:,2),'-k');



% Numerical Maps --------------------------------------------------------------
tau = 2*(1/12.5);
[tn,Sn] = delay_reconstruct(time,ankle_pos,(1/200),tau,3);

figure;
  plot(tn,Sn)

pause
numerical_poincare_ankleos(time,tn,Sn,locs,100);



return
index = 1:(numel(ankle_pos)-2*tau);

sn = ankle_pos(index);
sn_1tau = ankle_pos(index + tau);
sn_2tau = ankle_pos(index + 2*tau);

Sn = [sn, sn_1tau, sn_2tau];

figure;
for i=1:numel(sn)
  clf,
  plot3(sn(1:i), sn_1tau(1:i), sn_2tau(1:i),'-k'), hold on
  plot3(sn(i), sn_1tau(i), sn_2tau(i),'ok')
  drawnow
end


return

for i=1:numel(tau)
  fit{i} = numerical_poincare(time,ankle_pos,T,tau(i));
  alpha_vec = [alpha_vec, fit{i}.a];
end


figure, hold all
  plot(tau,alpha_vec,'-ok')
  plot(tau, ones(size(tau))*mean(alpha_vec),'--g')
  title('Numerical Agreement', ...
       'interpreter','latex', ...
       'fontsize', 24)
  xlabel('$\tau (s)$', ...
        'interpreter','latex', ...
        'fontsize', 20)
  ylabel('$a$', ...
        'interpreter','latex', ...
        'fontsize', 20)
  lh = legend('Numerical','Numerical Mean','Analytical');
  set(lh, ...
      'interpreter','latex', ...
      'fontsize', 15, ...
      'box', 'off', ...
      'Location','best')
grid on

figure, hold all
legend_str = {}
ColOrd = get(gca,'ColorOrder');
[m,n] = size(ColOrd);
colorIndex = 1;
for i=1:numel(tau)
  if (mod(i,10) == 0)
    ColRow = rem(colorIndex,m);
    if ColRow == 0
      ColRow = m;
    end
    Col = ColOrd(ColRow,:);
    plot(fit{i}.t_k,fit{i}.x_k,'o','Color',Col)
    plot(fit{i}.t_k(2:end),fit{i}.func(fit{i}.x_k(1:end-1)), '-','Color',Col)
    colorIndex = colorIndex + 1;
  end
end
  title('Fit: $x_{k+1} = x_k e^{aT} + x^* (1 - e^{aT})$', ...
        'interpreter','latex', ...
        'fontsize', 24)
 xlabel('Time (s)', ...
        'interpreter','latex', ...
        'fontsize', 20)
 ylabel('$x_k$', ...
        'interpreter','latex', ...
        'fontsize', 20)
 grid on

figure, hold all
legend_str = {}
ColOrd = get(gca,'ColorOrder');
[m,n] = size(ColOrd);
colorIndex = 1;
plot(time,ankle_pos,'-k')
for i=1:numel(tau)
  if (mod(i,10) == 0)
    ColRow = rem(colorIndex,m);
    if ColRow == 0
      ColRow = m;
    end
    Col = ColOrd(ColRow,:);
    plot(fit{i}.t_k,fit{i}.x_k,'o','Color',Col)
    colorIndex = colorIndex + 1;
  end
end
  title('Maps', ...
        'interpreter','latex', ...
        'fontsize', 24)
 xlabel('Time (s)', ...
        'interpreter','latex', ...
        'fontsize', 20)
 ylabel('$x(t)$', ...
        'interpreter','latex', ...
        'fontsize', 20)
 grid on


