% -----------------------------------------------------------------------------
% Flexi-force calibration data
%
% April 2015
%
% -----------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Sensor 1
% -----------------------------------------------------------------------------
S1.force_lbf = [0, 2.08, 12.10, 57.0, 113, 153, 195, 251];
S1.resistance_kOhms = [45000, 850, 125, 18, 7.9, 5.66, 4.46, 3.48];

% System Knowns
Rs = S1.resistance_kOhms*1000;
R = 270e3;
Vref = 1.65;
imax = 2.5e-3;

Vout = (R ./ (R + Rs)) * Vref;
iRs = (Vref - Vout) ./ (Rs);
iR = Vout ./ R;
i = (Vref) ./ (Rs + R);

figure; hold all;
  title('Resitive Response','fontsize',20);
  plot(S1.force_lbf, Rs, '*k')
  set(gca,'xscale','log','yscale','log');
  xlabel('Force (lbf)','fontsize',20);
  ylabel('Resitance (\Omega)','fontsize',20);
  grid on

figure; hold all;
  title('Voltage Response','fontsize',20);
  plot(S1.force_lbf, Vout, '-k')
  plot(S1.force_lbf, Vout.*0 + Vref,'--r')
  xlabel('Force (lbf)','fontsize',20);
  ylabel('Vout (V)','fontsize',20);
  grid on

figure; hold all;
  title('Current Response','fontsize',20);
  plot(S1.force_lbf, iRs, '-k')
  plot(S1.force_lbf, iR, '--r')
  plot(S1.force_lbf, i, '-.g')

  %plot(S1.force_lbf, i.*0 + imax,'--r')
  xlabel('Force (lbf)','fontsize',20);
  ylabel('I (A)','fontsize',20);
  grid on


return

figure; hold all;
  plot(S1.force_lbf, i_Rs);
  plot(S1.force_lbf, i_R,'--');
  plot(S1.force_lbf, i_max + i_R.*0,'--k');
  title('Current Response','fontsize',20);
  xlabel('Force (lbf)','fontsize', 20);
  ylabel('Current (mA)','fontsize',20);
  grid on

return

s = -R*Vref ./ (Rs + R).^2

figure;
  plot(S1.force_lbf, s )
  title('Sensitivity','fontsize',20);
  ylabel('V/lbs','fontsize',20);

return

R = 33;
Rs = S1.resistance_kOhms;
Vadc = 1.8;

% Simulate Voltage Divider 1
Vout1  = (R ./ (R + Rs)) * Vadc;
P_s1 = (Vadc - Vout1).^2 ./ (Rs*1000);
P_r1 = Vout1.^2 ./ (R*1000);




return

% Simulate Voltage Divider 2
Vout2  = (Rs ./ (R + Rs)) * Vadc;
P_s2 = Vout2.^2 ./ (Rs*1000);
P_r2 = (Vadc - Vout2).^2 ./ (R*1000);

figure;
  subplot(211); hold all;
    plot(S1.force_lbf, Vout1, '-k')
    plot(S1.force_lbf, Vout2,'--k')
    plot(S1.force_lbf, S1.force_lbf.*0 + Vadc,'--r')
    xlabel('Force (lbs)','fontsize',20)
    ylabel('Vout','fontsize',20)
    legend('Vout1', 'Vout2', 'Vadc')
    grid on
  subplot(212); hold all;
    plot(S1.force_lbf, P_s1,'-k');
    plot(S1.force_lbf, P_r1,'-r');
    plot(S1.force_lbf, P_s1 + P_r1,'-g');
    plot(S1.force_lbf, P_s2,'--ok');
    plot(S1.force_lbf, P_r2,'--or');
    plot(S1.force_lbf, P_s2 + P_r2,'--og');
    xlabel('Force (lbs)','fontsize',20)
    ylabel('Power (W)','fontsize',20)
    legend('Sensor Power 1', 'Resitor Power 1', 'Total Power 1', ...
           'Sensor Power 2', 'Resitor Power 2', 'Total Power 2', ...
           'location', 'SE')
    grid on

  figure; hold all;
    plot(S1.force_lbf, P_s1,'-k');
    plot(S1.force_lbf, P_r1,'-r');
    plot(S1.force_lbf, P_s1 + P_r1,'-g');
    figure; hold all;
    plot(S1.force_lbf, P_s2,'--k');
    plot(S1.force_lbf, P_r2,'--r');
    plot(S1.force_lbf, P_s2 + P_r2,'--g');
 
return
% -----------------------------------------------------------------------------
% Sensor 2
% -----------------------------------------------------------------------------
S2.force_lbf = [0, 10, 54, 100, 155, 207, 250];
S2.resistance_kOhms = [44000, 113, 16.58, 9.05, 5.82, 4.47, 3.79];

% -----------------------------------------------------------------------------
% Sensor 3
% -----------------------------------------------------------------------------
S3.force_lbf = [2, 12, 62, 103, 156, 215, 260];
S3.resistance_kOhms = [870, 290, 13.66, 7.86, 5.15, 3.72, 3.10];

% -----------------------------------------------------------------------------
% Sensor 4
% -----------------------------------------------------------------------------
S4.force_lbf = [2.5, 14.5, 50, 100, 165, 215, 263];
S4.resistance_kOhms = [730, 92.3, 21, 9.68, 5.97, 4.62, 3.9];

% -----------------------------------------------------------------------------
% Calibration Plot
% -----------------------------------------------------------------------------

figure; hold all;
  plot(S1.force_lbf,S1.resistance_kOhms,'ok')
  plot(S2.force_lbf,S2.resistance_kOhms,'or')
  plot(S3.force_lbf,S3.resistance_kOhms,'oc')
  plot(S4.force_lbf,S4.resistance_kOhms,'og')
  title('Force Sensor Calibration','interpreter','latex','fontsize',20);
  grid on
