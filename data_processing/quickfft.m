function rtn = quickfft(x,fs)
  L = numel(x);

  % Make L even
  if(mod(L,2) ~= 0)
    x = x(1:end-1);
    L = numel(x);
  end

  f = fs*(0:(L/2))/L;
  Y = fft(x)/L;
  Y = Y(1:(L/2)+1);

  figure;
    title('Spectrum','fontsize',20)
    plot(f,abs(Y),'k')
    xlim([0,0.1*fs])
    grid on

  figure; hold all;
    title('Reconstruction','fontsize',20)
    plot(x,'k')
    plot(ifft([Y; conj(flipud(Y(2:end-1)))]*L),'--r')


end
