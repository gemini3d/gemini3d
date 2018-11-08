function [thetat,phit]=geog2geomag(lat,lon)

%validateattr(lat, {'numeric'}, {'vector'}, mfilename)   %these can also be scalar
%validateattr(lon, {'numeric'}, {'vector'}, mfilename)   %can also be scalar
  
  thetan=11*pi/180;
  phin=289*pi/180;


  %enforce a [0,360] longitude
  inds=find(lon<0);
  loncorrected=lon;
  loncorrected(inds)=lon(inds)+360;
  inds=find(lon>360);
  loncorrected(inds)=lon(inds)-360;

  thetagp=pi/2-lat*pi/180;
%  phig=lon*pi/180;
  phig=loncorrected*pi/180;

  thetat=acos(cos(thetagp).*cos(thetan)+sin(thetagp).*sin(thetan).*cos(phig-phin));
  argtmp=(cos(thetagp)-cos(thetat).*cos(thetan))./(sin(thetat).*sin(thetan));
  alpha=acos( max(min(argtmp,1),-1) );
  phit=zeros(size(lat));
  inds=find( (phin>phig & phin-phig>pi) | (phin<phig & phig-phin<pi) );
  phit(inds)=pi-alpha(inds);
  inds=find( ~( (phin>phig & phin-phig>pi) | (phin<phig & phig-phin<pi) ) );
  phit(inds)=alpha(inds)+pi;

end
