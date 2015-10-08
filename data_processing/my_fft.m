function freqs = my_fft(x,dt,numOfFreq,doPlot)
  Fs = 1/dt;
  T = 1/Fs;
  L = numel(x);


  NFFT = 2^nextpow2(L);
  Y = fft(x(:,1),NFFT)/L;
  f = Fs/2*linspace(0,1,NFFT/2+1);
  MAG = 2*abs(Y(1:NFFT/2+1));

  [pks,loc] = findpeaks(MAG);

  [B,I] = sort(pks);
  pks = pks(I);
  loc = loc(I);

  if(doPlot)
  figure, hold all
    plot(f,MAG,'-k')
    for i=0:(numOfFreq-1)
      plot(f(loc(end-i*1)),MAG(loc(end-i*1)),'or')
    end
    title('Single-Sided Amplitude Spectrum of x(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
  end
    freqs = fliplr(f(loc(end-(numOfFreq-1):end)));
end
