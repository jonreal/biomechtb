function S_k = adaptiveILC(yd_k, y_k, scalar_gain, S_km1, varargin)
% adaptiveILC Iterative learning control with adaptive learning gains.
%             Currently works only in 1-D.
%
%   First iteration:
%   S_0 = adaptiveILC(yd_k, y_k, scalar_gain, 0, 'init', [gain,maxharmonic])
%
%
%   All others:
%   S_k = adaptiveILC(yd_k,y_k,scalar_gain,S_km1)
%
%                         yd_k  -  k_th desired output
%                         y_k   -  k_th output
%                         scalar_gain - inverse model (constant only!)
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
%   S_0 = adaptiveILC(yd_0, y_0, scalar_gain, 0, 'init', [gain, maxharmonic]);
%   S_1 = adaptiveILC(yd_1, y_1, scalar_gain, S_0);
%   S_2 = adaptiveILC(yd_2, y_2, scalar_gain, S_1);
%   .
%   .
%   .
%   S_k = adaptiveILC(yd_k, y_k, scalar_gain, S_km1);

  init = 0;
  nVarArgs = length(varargin);
  for i=1:2:nVarArgs
    switch varargin{i}
      case 'init'
        init = 1;
        gain = varargin{i+1}(1);
        maxharmonic = varargin{i+1}(2);
      otherwise
        fprintf('\n%s option not found!\n',varargin{i});
        return
    end
  end

  % Defaults
  alpha = 2;

  % Exclude last point (ie 100% sample == 0% sample)
  yd_k = yd_k(1:end-1);
  y_k = y_k(1:end-1);
  L = numel(yd_k);
  f = (0:L/2);
  f = f(:);

  % FFTs
  Yd_k = fft(yd_k);
  Yd_k = Yd_k(1:L/2+1);

  Y_k = fft(y_k);
  Y_k = Y_k(1:L/2+1);

  % Error
  E_k = Yd_k - Y_k;
  e_k_1 = norm(yd_k - y_k,1);
  e_k_2 = norm(yd_k - y_k,2);

  % Initialize S_0 if flag is set
  if (init)

    % Init rho
    rho_k = gain .* ones(numel(f),1);
    rho_k(maxharmonic+2:end) = 0;

    % First learning iteration
    U_kp1 = rho_k .* scalar_gain .* E_k;
    u_kp1 = ifft([U_kp1; conj(flipud(U_kp1(2:end-1)))]);

    % S_0 struct.
    S_0.f = f;
    S_0.Yd_k = Yd_k;
    S_0.yd_k = yd_k;
    S_0.Y_k = Y_k;
    S_0.y_k = y_k;
    S_0.rho_k = rho_k;
    S_0.E_bar_k = zeros(numel(f),1);
    S_0.U_bar_k = zeros(numel(f),1);
    S_0.E_k = E_k;
    S_0.U_kp1 = U_kp1;
    S_0.u_kp1 = u_kp1;
    S_0.e_k_1 = e_k_1;
    S_0.e_k_2 = e_k_2;

    S_k = S_0;
    return
  end

  % Unwrap input struct.
  rho_km1 = S_km1.rho_k;
  E_bar_km1 = S_km1.E_bar_k;
  U_bar_km1 = S_km1.U_bar_k;
  U_k = S_km1.U_kp1;
  E_km1 = S_km1.E_k;

  % Find at what frequencies the magnitude of error has increased/decreased
  E_incr = abs(E_k) > abs(E_km1);
  E_decr = ~(E_incr);

  % Updates
  rho_k = (E_incr .* rho_km1) / alpha + (E_decr .* rho_km1);
  E_bar_k = (E_incr .* E_bar_km1) + (E_decr .* E_k);
  U_bar_k = (E_incr .* U_bar_km1) + (E_decr .* U_k);

  U_kp1 = U_bar_k + rho_k .* scalar_gain .* E_bar_k;
  u_kp1 = ifft([U_kp1; conj(flipud(U_kp1(2:end-1)))]);

  % Pack S_kp1 struct.
  S_k.f = f;
  S_k.Yd_k = Yd_k;
  S_k.yd_k = yd_k;
  S_k.Y_k = Y_k;
  S_k.y_k = y_k;
  S_k.rho_k = rho_k;
  S_k.E_bar_k = E_bar_k;
  S_k.U_bar_k = U_bar_k;
  S_k.E_k = E_k;
  S_k.U_kp1 = U_kp1;
  S_k.u_kp1 = u_kp1;
  S_k.e_k_1 = e_k_1;
  S_k.e_k_2 = e_k_2;

end
