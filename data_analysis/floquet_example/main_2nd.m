
% Inputs
dt = 0.001;
wn_hz = 10;
zeta = 0.5;
A0 = 10;
w_hz = 1;
x0 = 0.1;
v0 = 0.1;

T = 1/w_hz;
t = 0:dt:(10*T);

n_s = 2;
n_f = 5;


A = [0, 1; -(2*pi*wn_hz)^2, -2*zeta*2*pi*wn_hz];
B = [0; 1];

u_step = 100*[0*t(1:n_s/dt)';
          ones(size(t((n_s/dt+1):(n_f/dt))))';
          0*t((n_f/dt+1):end)'];
u = @(t_) interp1(t,u_step,t_);


diffeq =@(t,x)  A*x + B*A0*(1 + u(t))*sin(2*pi*w_hz*t);


% Numerical -------------------------------------------------------------------
[t,x] = ode45(diffeq,t,[x0;v0]);

% Analytical  -----------------------------------------------------------------
wn = 2*pi*wn_hz;
w = 2*pi*w_hz;

lamda_bar = -zeta*wn;
w_bar = wn*sqrt(1 - zeta^2);
c4 = A0*( (wn^2 - w^2)/( (wn^2 - w^2)^2 + (2*zeta*wn*w)^2 ) );
c3 = -A0*( (2*zeta*wn*w)/( (wn^2 - w^2)^2 + (2*zeta*wn*w)^2 ) );
c1 = x0 - c3;
c2 = (v0 - lamda_bar*c1 - w*c4)/w_bar;

xx =@(t) exp(lamda_bar*t).*(c1*cos(w_bar*t) + c2*sin(w_bar*t)) ...
         + c3*cos(w*t) + c4*sin(w*t);
dxx =@(t) lamda_bar*exp(lamda_bar*t).*(c1*cos(w_bar*t) + c2*sin(w_bar*t)) ...
         + exp(lamda_bar*t).*(-c1*w_bar*sin(w_bar*t) + c2*w_bar*cos(w_bar*t)) ...
         -c3*w*sin(w*t) + c4*w*cos(w*t);


if(0)
% FFT -------------------------------------------------------------------------
freqs = my_fft(x,dt,2)

% Numerical Poincare  ---------------------------------------------------------
[tn, Sn] = delay_reconstruct(t,x(:,1),dt,T/5,2);

figure;
  plot(tn,Sn)

  return
alpha_vector = numerical_poincare_ankleos(t,tn,Sn,1:(1/1)/dt:numel(t),100)
end



% RLS -------------------------------------------------------------------------
phi_fun =@(t,w) [1, cos(t*w), sin(t*w), ...
                cos(2*t*w), sin(2*t*w), ...
                cos(3*t*w), sin(3*t*w), ...
                cos(4*t*w), sin(4*t*w), ...
                cos(5*t*w), sin(5*t*w), ...
                cos(6*t*w), sin(6*t*w), ...
                cos(7*t*w), sin(7*t*w), ...
                cos(8*t*w), sin(8*t*w)]';

window = 25;
lookAhead = T/8;

xd_meas = zeros(window,1);
t_meas = zeros(window,1);
xd_Hat = zeros(window,1);
theta0 = zeros(17,1);
P0 = 10*eye(17,17);
lam = 0.97;

theta_vector = zeros(numel(t),17);
yHat_vector = zeros(numel(t),1);

fh = figure;
  subplot(211),
     plot(0,0,'-k');
     plot(0,0,'--r');
  subplot(212)
     plot(0,zeros(1,17));

for i=1:numel(t)

  phi = phi_fun(t(i),1/T);
  [yHat, thetaHat, P] = rls(17,theta0,P0,lam,x(i,1),phi);

  theta0 = thetaHat;
  P0 = P;

  theta_vector(i,:) = theta0;
  yHat_vector(i) = yHat;

  if(0)
  figure(fh), clf
    subplot(211),
       plot(t(1:i),x(1:i),'-k');
       plot(t(1:i),yHat_vector(1:i),'--r');
    subplot(212)
       plot(t(1:i),theta_vector(1:i,:));
  drawnow;
  end
end

figure;
  plotyy(t,theta_vector',t,u(t))

figure;
  subplot(211)
    plot(t,x(:,1),'-k'), hold on
    plot(t,yHat_vector,'--r')
  subplot(212)
    plot(t,x(:,2),'-k'), hold on
    plot(t,dxx(t),'--r');

figure;
  plot(x(:,1),x(:,2),'-k');



