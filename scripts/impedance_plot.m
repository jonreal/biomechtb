if(~exist('c1_1'))
  c1_1 = vicon_process_trial('PA_C01_SS_F0_01');
end

force = c1_1.id.LAnkleMoment.X;
velocity = c1_1.id.LAnkleAngles.X;

figure;
  plot(velocity,force);

  pause

fs = 120;
L = numel(force);
f = fs*(0:(L/2))/L;

F = fft(force)/L;
V = fft(velocity)/L;
I = F.*V;

F = F(1:(L/2)+1);
V = V(1:(L/2)+1);
I = I(1:(L/2)+1);

figure; hold all;
  plot(f,abs(F),'k');
  plot(f,abs(V),'r');
  plot(f,abs(I),'g');
  title('Spectrum','fontsize',20);
  grid on

figure; hold all;
  plot(ifft([I; conj(flipud(I(2:end-1)))]*L),'--r')
  plot(force)
