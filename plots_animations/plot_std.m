function plot_std(t,X,c)
  % X is matrix in form [+std mean -std]

  t_fill = repmat([t'; flipud(t')], 1, 1);
  X_fill = repmat([X(:,1); flipud(X(:,3))], 1, 1);
  gca;
    hf = fill(t_fill, X_fill, c);
    set(hf,'EdgeColor','None','FaceAlpha',0.5)
    set(gca,'box','on')
end
