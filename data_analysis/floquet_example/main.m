% First order system dx/dt = ax + bsint

% Inputs ----------------------------------------------------------------------
a = -0.4;
b = 1;
w_hz = 1;
dt = 0.001;
x0 = 0;
snr = 45;

n_s = 2;
n_f = 2.5;



T = 1/w_hz;
t = 0:dt:10*T;

u_step = 0*[0*t(1:n_s/dt)';
          ones(size(t((n_s/dt+1):(n_f/dt))))';
          0*t((n_f/dt+1):end)'];
u = @(t_) interp1(t,u_step,t_);

% Numerical Sim ---------------------------------------------------------------
[t,x] = ode45(@(t,x) ...
              a*x + b*sin(2*pi*w_hz*t) + u(t), ...
              t, x0);

x_no_noise = x;
x = awgn(x,snr);


if(1)
% Delay Reconstruction --------------------------------------------------------
tau = 50;
index = 1:(numel(x)-tau);

sn = x(index);
sn_1tau = x(index + tau);

figure;
  plot(sn,sn_1tau)

figure;
  plot(x)
end

% Numerical Maps --------------------------------------------------------------
tau = linspace(0,T,50);
fit = {};
for i=1:numel(tau)
  fit{i} = numerical_poincare(t,x,T,tau(i));
end

% Analytical Poincare Maps ----------------------------------------------------
w = 2*pi*w_hz;

c1 = (x0 + b*w/(w^2 + a^2));
c2 = (-b*w/(w^2 + a^2));
c3 = (-a*b/(w^2 + a^2));

xx =@(t) c1*exp(a*t) + c2*cos(w*t) + c3*sin(w*t);

P_tau =@(x_k, tau, T) x_k*exp(a*T) ...
                   + (c2*cos(w*tau) + c3*sin(w*tau))*(1 - exp(a*T));
x_star =@(tau) (c2*cos(w*tau) + c3*sin(w*tau));


analytic = {};
for i=1:numel(tau)
  analytic{i}.t_k = tau:T:t(end);
  analytic{i}.x_k = interp1(t,x,analytic{i}.t_k,'cubic');
  analytic{i}.P_tau = P_tau(analytic{i}.x_k(1),tau(i),T);
  analytic{i}.x_star = x_star(tau(i));
end

% Compare Maps ----------------------------------------------------------------
alpha_vec = [];
for i=1:numel(tau)
  fprintf('--------------------------------------\n')
  fprintf('Tau = %f\n', tau(i))
  fprintf('Numerical\n')
  fprintf('a = %f, x* = %f\n', fit{i}.a, fit{i}.x_star)
  fprintf('Analytical\n')
  fprintf('a = %f, x* = %f\n', a, analytic{i}.x_star)
  fprintf('--------------------------------------\n\n')

  alpha_vec = [alpha_vec, fit{i}.a];
end

figure, hold all
  plot(t,x,'k')
  plot(t,x_no_noise,'--r')
  title('$\dot x = a x + b \sin \omega t$', ...
        'interpreter','latex', ...
        'fontsize', 24)
 xlabel('Time (s)', ...
        'interpreter','latex', ...
        'fontsize', 20)
 ylabel('$x(t)$', ...
        'interpreter','latex', ...
        'fontsize', 20)
  lh = legend('Noise','No Noise');
  set(lh, ...
      'interpreter','latex', ...
      'fontsize', 15, ...
      'box','off', ...
      'Location', 'best')
 grid on

e = (a - mean(alpha_vec))
figure, hold all
  plot(tau,alpha_vec,'-ok')
  plot(tau, ones(size(tau))*mean(alpha_vec),'--g')
  plot(tau, a*ones(size(tau)),'--r')
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
  set(gca,'ylim',[a-1,a+1])


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
plot(t,x,'-k')
for i=1:numel(tau)
  if (mod(i,5) == 0)
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
    plot(fit{i}.x_k(1:end-1) - fit{i}.x_star, ...
         fit{i}.x_k(2:end) - fit{i}.x_star, ...
         'o','Color',Col)
    plot(fit{i}.x_k(1:end-1) - fit{i}.x_star, ...
         fit{i}.func(fit{i}.x_k(1:end-1)) - fit{i}.x_star, ...
         '-','Color',Col)
    colorIndex = colorIndex + 1;
  end
end
  title('Fit: $x_{k+1} = x_k e^{aT} + x^* (1 - e^{aT})$', ...
        'interpreter','latex', ...
        'fontsize', 24)
 xlabel('$x_k - x^*$', ...
        'interpreter','latex', ...
        'fontsize', 20)
 ylabel('$x_{k+1} - x^*$', ...
        'interpreter','latex', ...
        'fontsize', 20)
 grid on

