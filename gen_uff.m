function [uff_k_plus_1, Ginv_k_plus_1] ...
           = gen_uff(uff_k, trial_k, bioIsRightFoot);

  fs = 120;

  if bioIsRightFoot
    yd = trial_k.Model.RAnkleMoment.X;
    bioHs = trial_k.Events.R.HS;

    y = trial_k.Model.LAnkleMoment.X;
    proHs = trial_k.Events.L.HS;

  else
    yd = trial_k.Model.LAnkleMoment.X;
    bioHs = trial_k.Events.L.HS;

    y = trial_k.Model.RAnkleMoment.X;
    proHs = trial_k.Events.R.HS;
  end

  fi = min(numel(bioHs),numel(proHs));
  yd = yd(bioHs(2):bioHs(fi));
  y = y(proHs(2):proHs(fi));

  if (mod(numel(yd),2) ~= 0)
    yd = yd(1:end-1);
  end
  if (mod(numel(y),2) ~= 0)
    y = y(1:end-1);
  end

  % Freq. Domain
  L = min(numel(yd), numel(y));
  f = (-L/2:L/2 - 1)*fs/L;

  Yd_k = fftshift(fft(yd,L));
  Y_k = fftshift(fft(y,L));
  Uff_k = fftshift(fft(uff_k,L));

  Ginv_k_plus_1 = (Uff_k./Y_k);
  Ginv_k_plus_1(isnan(Ginv_k_plus_1)) = 0;
  Ginv_k_plus_1(isinf(Ginv_k_plus_1)) = 0;

  Uff_k_plus_1 = Ginv_k_plus_1.*Y_k;

  figure;
    plot(f,abs(Uff_k_plus_1));

  uff_k_plus_1 = ifft(ifftshift(Uff_k_plus_1));

  figure; hold all;
    plot(uff_k_plus_1);
    plot(yd);
end
