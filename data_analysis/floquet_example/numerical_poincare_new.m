function alpha_vector = numerical_poincare_new(Sn,peakLocs,numMaps)

  debug = 0;

  % k-number of samples, n-states
  [k,n] = size(Sn);

  % Map locations indices
  peakLocs(peakLocs>k) = [];
  T_N = numel(peakLocs);

  % Create query points
  % Spit each period into numMaps points.
  for i=1:(T_N-1)
    iq(:,i) = linspace(peakLocs(i),peakLocs(i+1),numMaps);
  end
  alpha_vector = zeros(numMaps,n);
  A_tau = zeros(n,n,numMaps);
  G_tau = zeros(numMaps,n);
  maxFm_vector = zeros(numMaps,1);
  meanFM_vector = zeros(numMaps,1);
  stdFm_vector  = zeros(numMaps,1);

  fh = figure;
  fh1 = figure;
  for i=1:numMaps

    Sk = interp1(Sn,iq(i,:),'spline');

    % y-vector
    y = Sk(2:end,:);
    y = y(:);

    % A matrix
    kk = numel(Sk(1:end-1,1));
    G = kron(eye(n),ones(kk,1));
    TH = kron(eye(n),Sk(1:end-1,:));
    A = [TH, G];

    % Solve
    if (rank(A'*A) < kk)
      x_opt = pinv(A'*A)*A'*y;
    else
      x_opt = (A'*A)\A'*y;
    end
    x_opt(isnan(x_opt)) = 0;
    x_opt;

    % Eigenvalues
    A_tau(:,:,i) = reshape(x_opt(1:(n*n)),n,n);
    G_tau(i,:) = x_opt((n*n+1):end);
    FM = eig(A_tau(:,:,i));
    alpha_vector(i,:) = FM;

    maxFm_vector(i) = max(abs(FM));
    meanFm_vector(i) = mean(abs(FM));
    stdFm_vector(i) = std(abs(FM));

    if(debug)
      if (n >= 3)
      figure(fh);
        clf,
        subplot(121)
          plot3(Sn(:,1),Sn(:,2),Sn(:,3),'-k'), hold on
          plot3(Sk(:,1),Sk(:,2),Sk(:,3),'or')
          xlabel('$S_n$','interpreter','latex','fontsize',15)
          ylabel('$S_{n+\tau}$','interpreter','latex','fontsize',15)
          zlabel('$S_{n+2\tau}$','interpreter','latex','fontsize',15)
          title('Phase Space','interpreter','latex','fontsize',25)
          grid on
          box on
        subplot(122)
          plot(Sn(:,1),'-k'), hold on
          plot(iq(i,:),Sk(:,1),'or'),
        pause(0.01)

      elseif n == 2
      figure(fh);
        clf,
        plot(Sn(:,1),Sn(:,2),'-k'), hold on
        plot(Sk(:,1),Sk(:,2),'or')
        xlabel('$S_n$','interpreter','latex','fontsize',15)
        ylabel('$S_{n+\tau}$','interpreter','latex','fontsize',15)
        title('Phase Space','interpreter','latex','fontsize',25)
        grid on
        box on
        pause(0.01)
      end
    end
  end


  figure; hold on
    plot(1:numMaps,maxFm_vector,'-k')
    plot(1:numMaps,meanFm_vector,'--k')
    plot(1:numMaps,meanFm_vector + stdFm_vector','-r')
    plot(1:numMaps,meanFm_vector - stdFm_vector','-r')

    ylim([0,1])
    xlabel('\% Gait','interpreter','latex','fontsize',15)
    ylabel('max FM','interpreter','latex','fontsize',15)
    title('Maximum FM','interpreter','latex','fontsize',25)

  figure, hold all
    for i=1:n
      plot(1:numMaps,max(abs(alpha_vector(:,i))),'--o')
    end

   figure, hold all
    for i=1:n
      plot(real(alpha_vector(:,i)),imag(alpha_vector(:,i)),'ok')
    end
    plot(cos(linspace(0,2*pi)), sin(linspace(0,2*pi)))

end
