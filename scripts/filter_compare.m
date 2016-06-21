if(~exist('S'))
  S = embedded_process_data('filterTest');
end


fs = 1000;
quickfft(S.amp1s3,fs);



% Current Filter Design, FIR order 50
fc = 0.05;
wc_norm = 2*fc/fs;
b = fir1(50,wc_norm,'low')
[h,f] = freqz(b,1,500,fs);

[~,ii] = min(abs( mag2db(abs(h)) - (-3)));
approx_fc = f(ii)

mag1 = mag2db(abs(h));
phase1 = rad2deg(unwrap(angle(h)));
resp1 = filter(b,1,S.amp1s3);

% First order IIR, sampled
fc = 13;
s = tf('s');
sys = 2*pi*fc/(s + 2*pi*fc);
sys_d = c2d(sys,(1/1000),'tustin');
[num,den] = tfdata(sys_d);
b = num{1};
a = den{1};
[h,f] = freqz(b,a,500,1000);

[~,ii] = min(abs( mag2db(abs(h)) - (-3)));
approx_fc = f(ii)

mag2 = mag2db(abs(h));
phase2 = rad2deg(unwrap(angle(h)));
resp2 = filter(b,a,S.amp1s3);

% 2nd order butterworth
fc = 13;

[b,a] = butter(2,2*fc/fs,'low')
sys_d = tf(b,a,1000)

[b,a] = butter(2,2*fc/fs,'low','s')
sys_c = tf(b,a)

return

[h,f] = freqz(b,a,500,1000);

[~,ii] = min(abs( mag2db(abs(h)) - (-3)));
approx_fc = f(ii)

mag3 = mag2db(abs(h));
phase3 = rad2deg(unwrap(angle(h)));
resp3 = filter(b,a,S.amp1s3);

% Plots
figure; hold all
  subplot(211); hold all
    semilogx(f,mag1,'g')
    semilogx(f,mag2,'m')
    semilogx(f,mag3,'y')
    title('Mag','fontsize',20)
    semilogx([fc,fc],ylim,'--k')
    grid on
    set(gca,'xscale','log')
    xlim([0,500])
  subplot(212); hold all;
    semilogx(f,phase1,'g')
    semilogx(f,phase2,'m')
    semilogx(f,phase3,'y')
    semilogx([fc,fc,],ylim,'--k')
    title('Phase','fontsize',20)
    set(gca,'xscale','log')
    grid on
    legend('FIR','IIR1','IIR2')
    xlim([0,500])

figure; hold all
  title('Response','fontsize',20)
  plot(S.time,S.amp1s3,'k')
  plot(S.time,resp1,'g')
  plot(S.time,resp2,'m')
  plot(S.time,resp3,'y')
  xlabel('Time','fontsize',20);
  ylabel('bits','fontsize',20);
  legend('Raw','FIR','IIR1','IIR2')
  grid on


