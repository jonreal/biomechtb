k_s = 178000;
k_e = 180000;

ankle_pos = S.ankle_pos(k_s:k_e);
offset = 14;

[pks, locs] = findpeaks(ankle_pos)

figure, hold all
  subplot(211)
  plot(S.Kp(k_s:k_e), 'k')
  set(gca,'ylim',[0 200], ...
          'xlim',[0 (k_e - k_s)])
  xlabel('Sample','interpreter', ...
         'latex','fontsize',15)
  ylabel('Propotional Gain, $K_p$','interpreter','latex', ...
         'fontsize',15)
  box on
  grid on

  indx = find(pks>-0.02);
  locs = locs(indx);
  pks = pks(indx);

  for i=2:numel(locs)
    if (mod(i,2) == 0)
      ankle_pos(locs(i-1):locs(i)) = -ankle_pos(locs(i-1):locs(i));
    end
  end


  subplot(212), hold on
  plot(rad2deg(ankle_pos) + offset, 'k')

  [pks, locs] = findpeaks(rad2deg(ankle_pos) + offset)
  indx = find(pks>15);
  locs = locs(indx);
  pks = pks(indx);
  plot(locs,pks,'or')

  locs(2:end) - locs(1:end-1)

  meanSamp = mean(locs(2:end) - locs(1:end-1))

  meanFreq = meanSamp*(1/200)

  set(gca, ...
          'xlim',[0 (k_e - k_s)], ...
          'ylim',[0 25])
  xlabel('Sample','interpreter', ...
         'latex','fontsize',15)
  ylabel('Angle Angle (deg)','interpreter','latex', ...
         'fontsize',15)
  box on
  grid on

