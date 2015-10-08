function fit = numerical_poincare(t,x,T,tau)
% fit = numerical_poincare(t,x,T,tau)
% Function fits x_k+1 = x_k e^aT + b to 1-D periodic data
% using least squares || y - Ax ||^2 where,
% y = x_k+1
% A = [x_k  1]
% x = [alpha, beta]'
% alpha = exp(aT)
% beta = x_star(1 - exp(aT) = x_star(1 - alpha)

  t_k = tau:T:t(end)


  x_k = interp1(t,x,t_k,'cubic');

  y = x_k(2:end);
  y = y(:);

  xx = x_k(1:end-1);
  xx = xx(:);

  A = [xx, ones(size(xx))];

  % Weighting Matrix
  zeta = 1;
  omega = 0:(numel(xx)-1);
  w = exp((zeta - 1)*omega);
  W = diag(w);


  % Weighted Least Squares
  x_opt = (A'*W'*W*A)\(A'*W'*W*y);

  fit.t_k = t_k;
  fit.x_k = x_k;
  fit.a = log(x_opt(1))/T;
  fit.x_star = x_opt(2)/(1 - x_opt(1));
  fit.func =@(x_k) x_k*exp(fit.a*T) + fit.x_star*(1 - exp(fit.a*T));

end

