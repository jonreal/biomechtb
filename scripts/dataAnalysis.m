oldtrial1 = readData('test1_1_21');

figure; hold all;
  subplot(211), hold all;
    plot(oldtrial1.amp1s1,'k')
    plot(oldtrial1.amp1s2,'r')
    plot(oldtrial1.amp1s3,'g')
    legend('Toe','Mid','Heel')
    title('Old Trial 1 - Amp 1')
  subplot(212), hold all;
    plot(oldtrial1.amp2s1,'k')
    plot(oldtrial1.amp2s2,'r')
    plot(oldtrial1.amp2s3,'g')
    legend('Toe','Mid','Heel')
    title('Old Trial 1 - Amp 2')

oldtrial2 = readData('PA_A02_SS_F0_02_1_14');

figure; hold all;
  subplot(211), hold all;
    plot(oldtrial2.amp1s1,'k')
    plot(oldtrial2.amp1s2,'r')
    plot(oldtrial2.amp1s3,'g')
    legend('Toe','Mid','Heel')
    title('Old Trial 1 - Amp 1')
  subplot(212), hold all;
    plot(oldtrial2.amp2s1,'k')
    plot(oldtrial2.amp2s2,'r')
    plot(oldtrial2.amp2s3,'g')
    legend('Toe','Mid','Heel')
    title('Old Trial 1 - Amp 2')

trial1 = readData('PA_A02_SS_F0_01');

figure; hold all;
  subplot(211), hold all;
    plot(trial1.amp1s1,'k')
    plot(trial1.amp1s2,'r')
    plot(trial1.amp1s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 1 - Amp 1')
  subplot(212), hold all;
    plot(trial1.amp2s1,'k')
    plot(trial1.amp2s2,'r')
    plot(trial1.amp2s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 1 - Amp 2')

trial2 = readData('PA_A02_SS_F0_02');

figure; hold all;
  subplot(211), hold all;
    plot(trial2.amp1s1,'k')
    plot(trial2.amp1s2,'r')
    plot(trial2.amp1s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 2 - Amp 1')
  subplot(212), hold all;
    plot(trial2.amp2s1,'k')
    plot(trial2.amp2s2,'r')
    plot(trial2.amp2s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 2 - Amp 2')

trial3 = readData('PA_A02_SS_F1_02');

figure; hold all;
  subplot(211), hold all;
    plot(trial3.amp1s1,'k')
    plot(trial3.amp1s2,'r')
    plot(trial3.amp1s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 3 - Amp 1')
  subplot(212), hold all;
    plot(trial3.amp2s1,'k')
    plot(trial3.amp2s2,'r')
    plot(trial3.amp2s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 3 - Amp 2')

trial4 = readData('PA_A02_SS_F2_02');

figure; hold all;
  subplot(211), hold all;
    plot(trial4.amp1s1,'k')
    plot(trial4.amp1s2,'r')
    plot(trial4.amp1s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 4 - Amp 1')
  subplot(212), hold all;
    plot(trial4.amp2s1,'k')
    plot(trial4.amp2s2,'r')
    plot(trial4.amp2s3,'g')
    legend('Toe','Mid','Heel')
    title('Trial 4 - Amp 2')



