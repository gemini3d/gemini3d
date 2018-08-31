function ne=chapmana(z,nm,z0,H)
  zref=(z-z0)./H; 
  ne=nm*exp(0.5*(1-zref-exp(-zref)));
  
  ibad=find(ne<1);
  ne(ibad)=1;
end