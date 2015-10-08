function [yHat, thetaHat, P] = recurLS(n,theta0,P0,lam,y,phi)
% function [yHat, thetaHat, P] = recurLS(n,theta0,P0,lam,y,phi)
%
% Inputs:
%     n       -   filter order
%     theta0  -   old estimate
%     P0      -   old covariance of prediction error
%     lam     -   forgetting factor
%     y       -   data samples
%     phi     -   regressor vector

  K = P0*phi/(lam + phi'*P0*phi);
  P = ((eye(n) - K*phi')*P0)/lam;
  thetaHat = theta0 + K*(y - phi'*theta0);

  yHat = phi'*thetaHat;

return
