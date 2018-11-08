function [lat,lon]=geomag2geog(thetat,phit)

%validateattr(thetat, {'numeric'}, {'vector'}, mfilename)
%validateattr(phit, {'numeric'}, {'vector'}, mfilename)

  thetan=11*pi/180;
  phin=289*pi/180;


  %enforce phit = [0,2pi]
  inds=find(phit>2*pi);
  phitcorrected=phit;
  phitcorrected(inds)=phit(inds)-2*pi;
  inds=find(phit<0);
  phitcorrected(inds)=phit(inds)+2*pi;

%  thetag2p=acos(cos(thetat).*cos(thetan)-sin(thetat).*sin(thetan).*cos(phit));
  thetag2p=acos(cos(thetat).*cos(thetan)-sin(thetat).*sin(thetan).*cos(phitcorrected));

  beta=acos( (cos(thetat)-cos(thetag2p).*cos(thetan))./(sin(thetag2p).*sin(thetan)) );
%  phig2=zeros(size(phit));
  phig2=zeros(size(phitcorrected));
%  inds=find(phit>pi);
  inds=find(phitcorrected>pi);
  phig2(inds)=phin-beta(inds);
%  inds=find(phit<=pi);
  inds=find(phitcorrected<=pi);
  phig2(inds)=phin+beta(inds);
  inds=find(phig2<0);
  phig2(inds)=phig2(inds)+2*pi;
  inds=find(phig2>=2*pi);
  phig2(inds)=phig2(inds)-2*pi;

  thetag2=pi/2-thetag2p;
  lat=thetag2*180/pi;
  lon=phig2*180/pi;

end
