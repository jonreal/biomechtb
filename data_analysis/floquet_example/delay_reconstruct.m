function [Sn] = delay_reconstruct(sn,tau,m)
% function Sn = delay_reconstruct(sn,tau,m)
%
% Delay reconstruction for scalar measurements
% signal sn should have fixed sampling time (if not interp first)
%
% Input:
%   sn - vecotor of scalar measurements %TO DO: N-D data
%   tau - lag time (index delay)
%   m - embedding dimension
% Output:
%   Sn - m-dimensional state vector

  [n,k] = size(sn); % n-number of samples, k-dimension of scalar measurements
  Sn = zeros(n - m*tau,m);
  size(Sn)
  for i=1:m
    Sn(:,i) = sn([1:(n - m*tau)] + (i-1)*tau);
  end
end
