matFilePath = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/', ...
                'Data/PR_SS/PR_SS.mat'];

load(matFilePath);


sn = PR_SS.Model.RAnkleAngles.X;
t = PR_SS.time;
locs = PR_SS.Events.RHeelStrike_indx;
locs = locs(1:end-1);

figure, hold all
  plot(t,sn,'-k')
  plot(t(locs),sn(locs),'or')



tau = 10*(1/120);
[tn,Sn] = delay_reconstruct(t,sn,(1/120),tau,5);
numerical_poincare_ankleos(t,tn,Sn,locs,100);



figure;
  plot(PR_SS.Stats.gaitCycle,PR_SS.Stats.RMean_std.RAnkleAngles.X(:,2))
