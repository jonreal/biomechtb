function makeSineLookUpTab(filename)
  T = 2*pi;
  dt = T/1000;
  t = 0:dt:(T-dt);
  y = sin(t);

  figure;
    plot(t,y,'k')

  f = fopen(filename,'w');
  for i=1:numel(t)
    fprintf(f,'%f\n', y(i));
  end
  fclose(f);
end
