%if(~exist('S'))
%  S = embedded_process_data('~/research/data/PA_C01_Dec_1_2015/PA_C01_SS_F0_01');
%end

S = embedded_process_data_loc('~/gitrepo_master/embedded/anklebotOS/datalog/test1')

fs = 1000;
Ts = 1/fs;
fc = 10;

% Continous Filter
[b,a] = butter(2,fc*2*pi,'low','s');
sys_c = tf(b,a);

sys_d = c2d(sys_c,Ts,'tustin');
[num,den] = tfdata(sys_d);
num = num{1}
den = den{1}

%[b,a] = butter(2,2*fc/fs,'low');
%sys_d2 = tf(b,a,Ts);
%[num,den] = tfdata(sys_d2);
%num = num{1}
%den = den{1}


% Discrete Version (Bilinear Transform)
bb = a(2);
aa = a(3);

num = [Ts^2*aa, ...
       2*Ts^2*aa, ...
       Ts^2*aa]

den = [4 + bb*Ts + aa*Ts^2, ...
       2*aa*Ts^2 - 8, ...
       4 - bb*Ts + aa*Ts^2]
num = num./den(1)
den = den./den(1)

% Check

sys_d = c2d(sys_c,Ts,'tustin')
[num,den] = tfdata(sys_d);
num = num{1}
den = den{1}

sample = S.amp1.s3;
filt = sample.*0;
vel = sample.*0;
filt2 = sample.*0;

xx = [0; 0; 0];
yy = [0; 0; 0];

xxx = [0; 0];
yyy = [0; 0];


for i=1:numel(S.time)

  for k=3:-1:2
    xx(k) = xx(k-1);
    yy(k) = yy(k-1);
  end

  xx(1) = sample(i);

  % Filter
  filt(i) = num(1)*xx(1) + num(2)*xx(2) + num(3)*xx(3) ...
            - den(2)*yy(2) - den(3)*yy(3);
  yy(1) = filt(i);

  % Vel
  xxx(2) = xxx(1);
  yyy(2) = yyy(1);

  xxx(1) = filt(i);

  vel(i) = 2000*xxx(1) - 2000*xxx(2) - yyy(2);
  yyy(1) = vel(i);

end

rise_thrs = 500;
fall_thrs = -500;

figure; hold all;
  plot(sample);
  plot(filt)

figure;
  plot(vel)

  return



figure; hold all
  plot(filt./max(abs(filt)),'k');
  plot(vel./max(abs(vel)),'r')
