function fval = fobjective(x_k,T,a,b)

  v1 = x_k(2:end);
  v2 = x_k(1:end-1);
  v1 = v1(:);
  v2 = v2(:);

  fval = norm(v1 - v2*exp(a*T) - b + b*exp(a*T))^2;

end
