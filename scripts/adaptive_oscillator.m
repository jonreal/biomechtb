tt = PA_SS.time;
FF = PA_SS.Model.RAnkleMoment.X;

FF = PA_SS.Model.LAnkleMoment.X;


%dt = 0.001;
%tt = 0:dt:1;
%FF = sin(2*pi*30*tt) + 0.4*sin(2*pi*50*tt) + 0.5*sin(2*pi*3*tt);

figure;
  plot(tt,FF);


F =@(t) interp1(tt, FF, t);

figure;
  fd = quick_fft(FF,120,1)

return
K = 1;
mu = 0.001;
x0 = [0.01; 0; 1];

% Nonlinear oscilator - Hopf osicllator
fx =@(t,x,y,omega) (mu - sqrt(x.^2 + y.^2).^2).*x - omega.*y + K*F(t);
fy =@(t,x,y,omega) (mu - sqrt(x.^2 + y.^2).^2).*y + omega.*x;
fo =@(t,x,y) K*F(t).*(y./sqrt(x.^2 + y.^2));

FX =@(t,x) [fx(t,x(1),x(2),x(3)); ...
            fy(t,x(1),x(2),x(3));
            fo(t,x(1),x(2))];

[T,Y] = ode45(@(t,x) FX(t,x),tt,x0);

figure; hold all
  plot(T,Y)
  plot(T,fd(1)*ones(numel(T)));
  legend('x','y','w')


figure;
  plot(Y(:,1),Y(:,2))
