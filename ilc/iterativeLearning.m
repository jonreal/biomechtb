function rtn = iterativeLearning(trial_k, S_k, gain, bioIsRightFoot)

  % Unwrap input struct.
  U_k = S_k.U_k;
  U_k_m_1 = S_k.U_k_m_1;
  E_k_m_1 = S_k.E_k_m_1;
  rho_k_m_1 = S_k.rho_k_m_1;

  % Sample freq. (T_0/1000) - Fundamental Period / 1000 points
  fs = 1000;

  % # of FFT points
  L = 1000;

  % FFT freq. vector
  f = (0:L/2)*(fs/L);
  f = f(:);

  % Time vector;
  gaitCycle = trial_k.stats.gaitCycle(:);
  gaitCycle = gaitCycle(1:end-1);

  % Which foot is prosethetic?
  if bioIsRightFoot
    yd = trial_k.stats.r.RAnkleMoment.X(:,2);
    y = trial_k.stats.l.LAnkleMoment.X(:,2);
  else
    yd = trial_k.stats.l.LAnkleMoment.X(:,2);
    y = trial_k.stats.r.RAnkleMoment.X(:,2);
  end

  % Exclude last point (ie 100% gait sample == 0% gait sample)
  yd = yd(1:end-1);
  y = y(1:end-1);

  % FFT
  Y_k = fft(y);
  Y_k = Y_k(1:L/2 + 1);


  Yd_k = fft(yd);
  Yd_k = Yd_k(1:L/2 + 1);

  % Error
  E_k = Yd_k - Y_k;

  % rho update
  rho_k = zeros(1,L/2+1);
  for i=1:L/2+1
    if abs(E_k(i)) > abs(E_k_m_1(i))
      rho_k(i) = 0.5*rho_k_m_1(i);
      U_k(i) = U_k_m_1(i);
    else
      rho_k(i) = rho_k_m_1(i);
    end
  end
  rho_k = rho_k(:);
  U_k = U_k(:);

  % Learning update
  U_k_p_1 = U_k + gain.*rho_k.*E_k;

  u_k_p_1 = ifft([U_k_p_1; conj(flipud(U_k_p_1(2:end-1)))]);

  % Wrap output struct.
  rtn.f = f;
  rtn.U_k = U_k_p_1;
  rtn.U_k_m_1 = U_k;
  rtn.E_k_m_1 = E_k;
  rtn.rho_k_m_1 = rho_k;
  rtn.Yd_k_m_1 = Yd_k;
  rtn.Y_k_m_1 = Y_k;

  rtn.gaitCycle = gaitCycle;
  rtn.yd_k_m_1 = yd;
  rtn.y_k_m_1 = y;
  rtn.e_2 = norm(yd - y,2);
  rtn.e_1 = norm(yd - y,1);

  rtn.u_k = u_k_p_1;
end
