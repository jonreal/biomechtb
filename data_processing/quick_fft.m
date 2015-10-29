function [freqs, coeff] = quick_fft(x,fs,numOfFreq,doPlot)

  L = numel(x);
  Y = fftshift(fft(x,L))/L;
  Py = abs(Y);
  f = (fs*(-(L/2):(L/2-1))/L)';

  size(f)
  size(Py)

  % Remove zero freq. coef.
  PP = Py;
  PP(f == 0) = 0;

  % Find coef. at freq.
  [val,loc] = findpeaks(PP);
  numOfFreq = min(2*numOfFreq, numel(loc));
  [~,i] = sort(val,'descend');
  loc = sort(loc(i(1:numOfFreq)));

  if(doPlot)
    gcf, hold all;
      plot(f(f==0),10*log10(Py(f==0)),'or')
      plot(f(loc),10*log10(Py(loc)),'or');
      plot(f,10*log10(Py),'k');
      title('Amplitude Spectrum', ...
            'fontsize',20)
      ylabel('10log_{10}|Y(f)|','fontsize',20)
      xlabel('Frequency (Hz)','fontsize',20);
      xlim([-f(loc(numOfFreq)) - 50*(f(2)-f(1)), ...
            f(loc(numOfFreq)) + 50*(f(2)-f(1))])
      grid on
  end
  freqs = f(loc(1:numOfFreq));
  coeff = Y(loc(1:numOfFreq));

  freqs = [0; freqs(:)];
  coeff = [Y(f==0); coeff(:)];
end
