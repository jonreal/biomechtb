function si = symmetryIndex(f1, f2)
  % SYMMERTYINDEX computes the symmetry index of two gait variables (e.g.,
  % left/right, prosthetic/intact). The gait variables should be of the form
  % f1 = f(%gait), e.g., 0-100% gait.
  %
  % si = symmetryIndex(f1,f2)
  %
  %   si - symmetry index
  %   f1 - gait variable 1
  %   f2 - gait variable 2

  si = (f1 - f2)./(0.5*(f1 + f2))*100;
end
