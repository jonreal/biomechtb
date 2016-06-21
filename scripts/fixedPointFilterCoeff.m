samplingFrequency = 1000;
centerFrequency = 100;
transitionWidth = 50;
passbandRipple = 1;
stopbandAttenuation = 80;

designSpec = fdesign.lowpass('Fp,Fst,Ap,Ast',...
      centerFrequency-transitionWidth/2, ...
      centerFrequency+transitionWidth/2, ...
      passbandRipple,stopbandAttenuation, ...
      samplingFrequency);
LPF = design(designSpec,'equiripple','SystemObject',true)

fvtool(LPF,'Analysis','phase')
coeff = round(LPF.Numerator * 32768);

numOfCoeff = numel(coeff);

fileName = input('\nEnter file name for filter coeff\n','s');

file = fopen(['./',fileName],'w');
for i=1:numel(coeff)
  fprintf(file,'%i\n', coeff(i));
end
fclose(file);

