matFilePath = ['~/Dropbox/Professional/UW_PHD/', ...
               'Prosthetic_Research/Data/PA_A01/', ...
                'Data/PR_SS/PR_SS.mat'];

load(matFilePath);

% Input data
%x = PR_SS.Markers.CentreOfMass.Z;
x = PR_SS.Model.RAnkleAngles.X;
t = PR_SS.time;

% Delay Reconstuction ---------------------------------------------------------
figure;
  autocorr(x,100);

x_corr = autocorr(x,100);
tau = find(x_corr<0);
tau = tau(1)

m = 3;
Sn = delay_reconstruct(x,tau,m);


Sn = [PR_SS.Model.RAnkleAngles.X, ...
      PR_SS.Model.d_RAnkleAngles.X, ...
      PR_SS.Model.d2_RAnkleAngles.X];
figure, hold all
  plot(Sn(:,1),'-k')
  title('Data Set','interpreter','latex','fontsize',25)
  xlabel('Sample \#','interpreter','latex','fontsize',15)
  ylabel('Ankle Angle (rad)','interpreter','latex','fontsize',15)

figure, hold all
  plot3(Sn(:,1),Sn(:,2),Sn(:,3),'-k')
  title('Data Set','interpreter','latex','fontsize',25)
  xlabel('Sample \#','interpreter','latex','fontsize',15)
  ylabel('Ankle Angle (rad)','interpreter','latex','fontsize',15)
  view(3)
% --- FFT ---------------------------------------------------------------------
w = my_fft(PR_SS.Model.RAnkleAngles.X,1,1,0);
meanT = ceil(1/w);

% --- Poincare ----------------------------------------------------------------
if(1)
  numOfMaps = 1001;
  heelStrike = PR_SS.Events.LHeelStrike_indx;
  heelStrike(heelStrike>numel(Sn(:,1))) = [];
  [maxFM,FM] = floquet_multipliers(t,Sn,numOfMaps,heelStrike);

  figure, hold on
    plot(linspace(0,100,numOfMaps),maxFM,'-ok')
    plot(linspace(0,100,numOfMaps),mean(maxFM).*ones(numOfMaps,1))
    xlim([0,100])
    ylim([0,1])
end

% -- Lyapunov Exponent --------------------------------------------------------
if(1)
  d = lyapunov_exponent(Sn,meanT,4,1);

end


return

