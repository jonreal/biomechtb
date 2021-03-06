function S_k = adaptiveILC(yd_k, y_k, gamma_k, gain_k, S_km1, varargin)
% adaptiveILC Iterative learning control with adaptive learning gains.
%             Currently works only in 1-D.
%
%
%   Algorithm:
%
%     at each harmonic (in freq domain)
%
%       condition: |E_k| > ( |E_bar_km1| + gamma_k )
%
%                        / rho_km1 alpha  if condition
%               rho_k = {
%                        \ rho_km1
%
%                          / E_bar_km1  if condition
%               E_bar_k = {
%                          \ E_k
%
%                          / U_bar_km1  if condition
%               U_bar_k = {
%                          \ U_k
%
%               U_kp1 = U_bar_k + rho_k * gain_k * E_bar_k
%
%
%   Note: rho_k always is in a range [0,1], use gain_k to set the learning gain
%   gamma_k is padding for error comparison. Use std of output signal druing
%   k=0 iteration.
%
%   Use:
%
%   First iteration:
%   S_0 = adaptiveILC(yd_k, y_k, gain_k, 0, 'init', maxharmonic)
%
%
%   All others:
%   S_k = adaptiveILC(yd_k, y_k, gain_k, S_km1)
%
%                         yd_k    -  k_th desired output
%                         y_k     -  k_th output
%                         gamma_k -  error buffer (use variance of output k=0)
%                         gain_k  -  learning weight
%                         S_km1   -  {k-1}_th stucture
%
%                         S_k : Yd_k     - k_th desired output
%                               yd_k     - k_th desired output (time domain)
%                               Y_k      - k_th ouput
%                               y_k      - k_th ouput (time domain)
%                               rho_k    - k_th learning gains
%                               E_bar_k  - k_th update error
%                               U_bar_k  - k_th update input
%                               E_k      - k_th error
%                               U_kp1    - (k+1)_th input
%                               u_kp1    - (k+1)_th input (time domain)
%                               e_k_1    - k_th norm(yd_k - y_k,1)
%                               e_k_2    - k_th norm(yd_k - y_k,2)
%
% Example:
%   S_0 = adaptiveILC(yd_0, y_0, gain_k, 0, 'init', maxharmonic);
%   S_1 = adaptiveILC(yd_1, y_1, gain_k, S_0);
%   S_2 = adaptiveILC(yd_2, y_2, gain_k, S_1);
%   .
%   .
%   .
%   S_k = adaptiveILC(yd_k, y_k, gain_k, S_km1);

  init = 0;
  nVarArgs = length(varargin);
  for i=1:2:nVarArgs
    switch varargin{i}
      case 'init'
        init = 1;
        maxharmonic = varargin{i+1};
      otherwise
        fprintf('\n%s option not found!\n',varargin{i});
        return
    end
  end

  % Defaults
  alpha = 2;
  div_const = 1;

  % Number of points
  L = numel(yd_k);

  % This is really f/fs -> fs is the sampling freq.
  % Typically use normalized data (e.g., %gait) hence fs = NNpoints/T_gait
  f = (0:(L/2))/L;
  f = f(:);

  % FFTs
  Yd_k = fft(yd_k);
  Yd_k = Yd_k(1:(L/2+1));

  Y_k = fft(y_k);
  Y_k = Y_k(1:(L/2+1));

  % Error
  E_k = Yd_k - Y_k;

  e_k_inf = norm(yd_k - y_k,inf);
  e_k_2 = norm(yd_k - y_k,2);
  E_k_inf = norm(Yd_k - Y_k,inf);
  E_k_2 = norm(Yd_k - Y_k,2);

  % Initialize S_0 if flag is set
  if (init)

    % Init rho
    rho_k = ones(numel(f),1);
    rho_k((maxharmonic+2):end) = 0;

    U_kp1 = rho_k .* gain_k .* E_k;
    u_kp1 = ifft([U_kp1; conj(flipud(U_kp1(2:end-1)))]);

    % S_0 struct.
    S_0.f = f;
    S_0.Yd_k = Yd_k;
    S_0.yd_k = yd_k;
    S_0.Y_k = Y_k;
    S_0.y_k = y_k;
    S_0.rho_k = rho_k;
    S_0.gain_k = gain_k;
    S_0.gamma_k = gamma_k;

    S_0.E_bar_k = E_k;
    S_0.U_bar_k = zeros(numel(f),1);

    % To monitor growth of reference signal
    S_0.E_0 = E_k;
    S_0.Yd_0 = Yd_k;
    S_0.update = E_k.*0;

    S_0.E_k = E_k;
    S_0.U_kp1 = U_kp1;
    S_0.u_kp1 = u_kp1;

    S_0.e_k_inf = e_k_inf;
    S_0.e_k_2 = e_k_2;
    S_0.E_k_inf = E_k_inf;
    S_0.E_k_2 = E_k_2;

    S_k = S_0;
    return
  end

  % Unwrap input struct.
  Yd_km1 = S_km1.Yd_k;
  rho_km1 = S_km1.rho_k;
  U_k = S_km1.U_kp1;
  E_km1 = S_km1.E_k;

  E_bar_km1 = S_km1.E_bar_k;
  U_bar_km1 = S_km1.U_bar_k;

  % Check reference signal growth
  divergence = abs(S_km1.Yd_0 - Yd_k) > div_const*abs(S_km1.E_0);

  % Find at what frequencies the magnitude of error has increased/decreased
  E_incr = abs(E_k) > (abs(E_bar_km1) + gamma_k);
  E_decr = ~E_incr;

  % Update only if ~divergence and E_decr
  update = (~divergence) & E_decr;

  % Updates
  rho_k = ((~update).* rho_km1/alpha) + (update.* rho_km1);
  E_bar_k = ((~update).* E_bar_km1) + (update.* E_k);
  U_bar_k = ((~update).* U_bar_km1) + (update.* U_k);

  U_kp1 = U_bar_k + rho_k .* gain_k .* E_bar_k;
  u_kp1 = ifft([U_kp1; conj(flipud(U_kp1(2:(end-1))))]);

  % Copy previous structure
  S_k = S_km1;

  % Update
  S_k.update = update;
  S_k.Yd_k = Yd_k;
  S_k.yd_k = yd_k;
  S_k.Y_k = Y_k;
  S_k.y_k = y_k;
  S_k.E_bar_k = E_bar_k;
  S_k.U_bar_k = U_bar_k;
  S_k.E_k = E_k;
  S_k.U_kp1 = U_kp1;
  S_k.u_kp1 = u_kp1;
  S_k.e_k_inf = e_k_inf;
  S_k.e_k_2 = e_k_2;
  S_k.E_k_inf = E_k_inf;
  S_k.E_k_2 = E_k_2;
end
