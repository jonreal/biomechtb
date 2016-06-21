% Current Configuration - (Jan 13 2016)
ell = 0.004;
eff = 0.9;
R = 27/5;
k_t = 0.00553;

a0 = 0.1;
d = 0.06;
theta0 = deg2rad(-4.5);
b = sqrt(a0^2 + d^2);
beta0 = atan(a0/d);

beta_func =@(theta) beta0 + theta - theta0;
a_func =@(theta) sqrt(d^2 + b^2 - 2*d*b*cos(beta_func(theta)));
phi_func =@(theta) asin(d*sin(beta_func(theta))./a_func(theta)) ...
                    + beta_func(theta) - pi/2;

ankleTorque2current =@(theta) (1/k_t) ...
                     * (1/R) ...
                     * (ell/(2*pi*eff)) ...
                     * (1./(d*cos(phi_func(theta))));

theta_range = deg2rad(-20:20);
figure, hold all
  plot(rad2deg(theta_range), ankleTorque2current(theta_range))
  plot(rad2deg(theta_range), ...
       mean(ankleTorque2current(theta_range)) + 0.*theta_range)
  grid on

% Use 0.4 as scaling
mean(ankleTorque2current(theta_range))

