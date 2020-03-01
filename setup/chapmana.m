function ne=chapmana(z,nm,z0,H)
  zref=(z-z0)./H;
  ne=nm*exp(0.5*(1-zref-exp(-zref)));

  ne(ne<1) = 1;
end