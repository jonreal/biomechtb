function rtn = nthHarm(F, n)
  FF = F.*0;
  FF(n) = F(n);

  rtn.F_n = F(n);
  rtn.f_n = ifft([FF; conj(flipud(FF(2:end-1)))]);
end
