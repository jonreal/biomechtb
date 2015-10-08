% Params
sig = 16;
R = 45.92;
b = 4;
x0 =[-10,20,40]
x0 = [0.01,0.01,0.01]
x0 = [-14.1865  -13.0232   47.3509]

% Diff eq
f =@(t,x) [sig*(x(2) - x(1)); x(1)*(R-x(3)) - x(2); x(1)*x(2) - b*x(3)];

[t,x] = ode45(f,[0:0.01:(5000*0.01)],x0)



figure;
  subplot(311)
    plot(t,x(:,1),'-k')
    ylabel('$x_1$','interpreter','latex','fontsize',15)
    title('Lorenz Attractor - Time Response', ...
          'interpreter','latex','fontsize',15)
  subplot(312)
    plot(t,x(:,2),'-k')
    ylabel('$x_2$','interpreter','latex','fontsize',15)
  subplot(313)
    plot(t,x(:,3),'-k')
    ylabel('$x_3$','interpreter','latex','fontsize',15)
    xlabel('$t$ $(s)$','interpreter','latex','fontsize',15)

  figure;
    plot3(x(:,1),x(:,2),x(:,3),'k'),
    grid on
    title('Lorenz Attractor - Phase Space','interpreter','latex','fontsize',25)
    xlabel('$x_1$','interpreter','latex','fontsize',15)
    ylabel('$x_2$','interpreter','latex','fontsize',15)
    zlabel('$x_3$','interpreter','latex','fontsize',15)

sn = x(:,1);
% FFT
w = my_fft(sn,1,1,1)
meanT = ceil(1/w)

tau = 50;
m = 3;

if(1)
  figure;
    autocorr(sn,100)
  figure;
    [FNN] = fnn_deneme(sn,tau,10,15,2)
end


Sn = delay_reconstruct(sn,tau,m);

figure;
  plot3(Sn(:,1),Sn(:,2),Sn(:,3),'k')
  grid on
  title('Phase Space Reconstuction', ...
        'interpreter','latex','fontsize',25)
  xlabel('$x_1(t)$','interpreter','latex','fontsize',15)
  ylabel('$x_1(t + \tau)$','interpreter','latex','fontsize',15)
  zlabel('$x_2(t + 2\tau)$','interpreter','latex','fontsize',15)


