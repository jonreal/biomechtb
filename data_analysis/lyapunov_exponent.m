function [avg_dj,lyp_exp] = ...
            lyapunov_exponent(Sn,T,numOfStrides,strides2fit)
  debug=1;
  if(debug)
    fh = figure;
    fh2 = figure;
  end

  M = numel(Sn(:,1));     % # of data points
  I = numOfStrides*T;     % Number of strides (analysis time)

  j_hat = zeros(1,M);
  d0_j = zeros(1,M);

  % First find all NN
  for j=1:M

    % Constrain NN search
    qindx = [[1:(j-T)],[(j+T):M]];
    [indx,dst] = knnsearch(Sn(qindx,:),Sn(j,:));

    % Store nearest neighbor and distance
    j_hat(j)= qindx(indx);
    d_j(j) = norm(Sn(j,:) - Sn(qindx(indx),:),2);

    if((debug) && (mod(j,100) == 0))
      figure(fh), clf, hold all
        plot(Sn(:,1),'-k')
        plot(qindx,Sn(qindx),'--r')
        plot(j,Sn(j,1),'or')
        plot(j_hat(j),Sn(j_hat(j),1),'oc')
        pause(0.001)
    end

  end

  avg_dj = zeros(1,I);
  % For each NN track distant as a function of time (i)
  avg_dj(1) = mean(log(d_j));
  for i=1:I
    d_j = [];
    for j=1:M
      if ( ((j_hat(j)+i) <= M) && ((j+i) <= M) )
        d_j(end+1) = norm(Sn(j+i,:) - Sn(j_hat(j)+i,:),2);

        if(debug)
          figure(fh), clf, hold all
            plot(Sn(:,1),'-k')
            plot(j,Sn(j,1),'or')
            plot(j_hat(j),Sn(j_hat(j),1),'oc')
            plot(j+i,Sn(j+i,1),'*r')
            plot((j_hat(j)+i),Sn(j_hat(j)+i,1),'*c')
            pause(0.001)
        end
      end
    end
    avg_dj(i+1) = mean(log(d_j));
    numel_avg_dj(i+1) = numel(d_j);
  end

  xx = ((0:strides2fit*T)./T)';
  yy = avg_dj(1:numel(xx))';
  lyp_fit = fit(xx,yy,'poly1')
  lyp_exp = lyp_fit.p1;

  figure, hold all
    plot((0:I)./T,avg_dj,'-k')
    plot(lyp_fit,xx,yy,'--r')
    xlabel('\# strides','interpreter','latex','fontsize',15)
    ylabel('$<$ln$(d_j)>$','interpreter','latex','fontsize',15)
    title('Average Logarithmic Divergence','interpreter','latex','fontsize',25)
    grid on
    box on
    xlim([0,numOfStrides])
    yl = ylim;
    plot(strides2fit*ones(1,10),linspace(yl(1),yl(2),10),'--k')

   figure;
    plot(numel_avg_dj)
end

