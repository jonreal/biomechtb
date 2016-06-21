% Paramaters
I = 1;
b = 0.5;

kp = [41, 301, 617, 439]';
Tp =@(x) -sum(kp(1) + kp(2).*x + kp(3).*x.^2 + kp(4).*x.^3, 1);
Ta =@(t,x) 0;


f =@(t,x) [x(2); ...
         Tp(x(1)) - b*x(2)];
g =@(t,x) [0;
         Ta(t,x)];

xd =@(t,x) f(t,x) + g(t,x); 


tt = linspace(0,10,1000);

[T,X] = ode45(@(t,x) ...
              xd(t,x), tt, [0;0]);

y = Tp(X(:,1)) - b*(X(:,2)) + Ta(T,X(:,1));

figure; hold all
  plot(tt,y)
  plot(tt,Y(:,1))
