fs = 120;
tt = PA_SS.time;
gait = PA_SS.Stats.gaitCycle;

Tp = PA_SS.Model.LAnkleMoment.X;
Tb = PA_SS.Model.RAnkleMoment.X;

xp = PA_SS.Model.LAnkleAngles.X;
xb = PA_SS.Model.RAnkleAngles.X;

Pp = PA_SS.Model.LAnklePower.X;
Pb = PA_SS.Model.RAnklePower.X;

xp_stat = PA_SS.Stats.LMean_std.LAnkleAngles.X;
xb_stat = PA_SS.Stats.RMean_std.RAnkleAngles.X;

Tp_stat = PA_SS.Stats.LMean_std.LAnkleMoment.X;
Tb_stat = PA_SS.Stats.RMean_std.RAnkleMoment.X;

Pp_stat = PA_SS.Stats.LMean_std.LAnklePower.X;
Pb_stat = PA_SS.Stats.RMean_std.RAnklePower.X;

% Freq. Domain Analysis
L_hs = PA_SS.Events.L.HS;
R_hs = PA_SS.Events.R.HS;

fi = min(numel(L_hs),numel(R_hs));

Tp_aligned = Tp(L_hs(2):L_hs(fi));
Tb_aligned = Tb(R_hs(2):R_hs(fi));

if (mod(numel(Tp_aligned),2) ~= 0)
  Tp_aligned = Tp_aligned(1:end-1);
end
if (mod(numel(Tb_aligned),2) ~= 0)
  Tb_aligned = Tb_aligned(1:end-1);
end

figure; hold all;
  title('Ankle Torque','fontsize',20);
  plot(Tp_aligned,'k');
  plot(Tb_aligned,'r');
  xlabel('Sample','fontsize',20);
  ylabel('Torque (Nm/kg)','fontsize',20);
  legend('Pros','Bio');

% Freq. vectors
Lp = numel(Tp_aligned);
Lb = numel(Tb_aligned);

fp = (-Lp/2:Lp/2 - 1)*fs/Lp;
fp = fp(:);

fb = (-Lb/2:Lb/2 - 1)*fs/Lb;
fb = fb(:);

TTp = fftshift(fft(Tp_aligned));
TTb = fftshift(fft(Tb_aligned));

figure; hold all;
  stem(fp(fp >= 0),abs(TTp(fp >=0)),'k');
  stem(fb(fb >=0),abs(TTb(fb >=0)),'--r');
  xlim([0,10]);












