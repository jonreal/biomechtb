function S = embedded_process_data(file)
  D = dlmread(['./',file],'\t',6,0);

  S.dt = 1/1000;

  D(D(:,2) ~= 1,:) = [];

  S.timeStamp = D(:,1);
  S.time = (D(:,1) - D(1,1)).*1/1000;

  S.sync = D(:,2);
  S.avgT_cnts = D(:,3);
  S.heelStrikeCnt = D(:,4);
  S.gaitPhase = D(:,5);
  S.anklePos = D(:,6);
  S.ankleVel = D(:,7);
  S.fbCurrentCmd = D(:,8);
  S.ffCurrentCmd = D(:,9);
  S.motorDuty = D(:,10);
  S.motorCurrent = D(:,11);
  S.motorVel = D(:,12);
  S.amp1s1 = D(:,13);
  S.amp1s2 = D(:,14);
  S.amp1s3 = D(:,15);
  S.amp2s1 = D(:,16);
  S.amp2s2 = D(:,17);
  S.amp2s3 = D(:,18);

  S.acc1 = D(:,19);
  S.acc2 = D(:,20);
  S.acc3 = D(:,21);
  S.gyro1 = D(:,22);
  S.gyro2 = D(:,23);
  S.gyro3 = D(:,24);

end
