function [ns,Ts,vsx1] = eqICs3D(p, xg)
%% generate (arbitrary) initial conditions for a grid.
% NOTE: only works on symmmetric closed grids!
%
% [f107a, f107, ap] = activ;
%
narginchk(2,2)
validateattributes(p, {'struct'},{'scalar'})
validateattributes(xg, {'struct'},{'scalar'})

%% MAKE UP SOME INITIAL CONDITIONS FOR FORTRAN CODE
mindens=1e-100;

%% SLICE THE FIELD IN HALF IF WE ARE CLOSED
natm = msis_matlab3D(p, xg);

closeddip = abs(xg.r(1,1,1) - xg.r(xg.lx(1),1,1)) < 50e3;     %logical flag marking the grid as closed dipole
if closeddip         %closed dipole grid
%    [~,ialtmax]=max(xg.alt(:,1,1));
%    lalt=ialtmax;
  lalt=floor(xg.lx(1)/2);                         %FIXME:  needs to work with asymmetric grid...
  alt=xg.alt(1:lalt,:,:);
  lx1=lalt;
  lx2=xg.lx(2);
  lx3=xg.lx(3);
  Tn=natm(1:lalt,:,:,4);
  g=abs(xg.gx1(1:lalt,:,:));
  g=max(g,1);
  for ix3=1:lx3
    for ix2=1:lx2
      [~,ialt]=min(abs(g(:,ix2,ix3)-1));
      if ialt~=lalt
        g(ialt:lalt,ix2,ix3)=1;
      end
    end
  end
else
  alt=xg.alt;
  lx1=xg.lx(1);
  lx2=xg.lx(2);
  lx3=xg.lx(3);
  Tn=natm(:,:,:,4);
  g=abs(xg.gx1);
end


%CONSTANTS
kb=1.38e-23;
amu=1.67e-27;


ns=zeros(lx1,lx2,lx3,7);
for ix3=1:lx3
  for ix2=1:lx2
    Hf = kb * Tn(:,ix2,ix3) / amu /16 ./ g(:,ix2,ix3);
    z0f=325e3;
    He = 2 * kb * Tn(:,ix2,ix3) / amu / 30 ./ g(:,ix2,ix3);
    z0e = 120e3;
    ne = chapmana(alt(:,ix2,ix3), p.nmf,z0f,Hf)+chapmana(alt(:,ix2,ix3), p.nme,z0e,He);
    rho = 1/2*tanh((alt(:,ix2,ix3)-200e3)/45e3)-1/2*tanh((alt(:,ix2,ix3)-1000e3)/200e3);

    inds=find(alt(:,ix2,ix3)>z0f);
    if ~isempty(inds)
      n0 = p.nmf;
      %     [n0,ix1]=max(ne);  %in case it isn't exactly z0f
      %     if xg.r(1,1)>xg.r(2,1)
      %         inds=1:ix1;
      %     else
      %         inds=ix1:lx1;
      %     end
      ms = rho(inds)*16*amu+(1-rho(inds))*amu;      %topside composition only
      H=kb*2*Tn(inds,ix2,ix3)./ms./g(inds,ix2,ix3);
      z=alt(inds,ix2,ix3);
      lz=numel(z);
      [z,iord]=sort(z);
      %     z=[z; 2*z(lz)-z(lz-1)];
      z = [z0f; z];
      integrand = [1./H(iord)];
      integrand = [integrand; integrand(lz)];
%        redheight=intrap(integrand,z);
      redheight=cumtrapz(z,integrand);
      netop=n0*exp(-redheight);
      nesort=zeros(lz,1);
      for iz=1:lz
        nesort(iord(iz))=netop(iz);
      end
      ne(inds) = nesort;
    end


    %% O+
    ns(:,ix2,ix3,1)= rho .* ne;
    zref=900e3;
    inds0=find(alt(:,ix2,ix3)>zref);
    if ~isempty(inds0)
      [altsort,iord]=sort(alt(:,ix2,ix3));
      nsort=ns(:,ix2,ix3,1);
      nsort=nsort(iord);
%        n0=interpolate(nsort,altsort,zref,'lin','lin');
      n0=interp1(altsort,nsort,zref);
      %     [tmp,iref]=min(abs(alt(:,ix2,ix3)-900e3));
      %     if xg.r(1,1)>xg.r(2,1)
      %         inds0=1:iref;
      %     else
      %         inds0=iref:lx1;
      %     end
      %    n0=ns(iref,ix2,ix3,1);
      ms=16*amu;
      H=kb*2*Tn(inds,ix2,ix3)./ms./g(inds,ix2,ix3);
      z = alt(inds0,ix2,ix3);
      lz=numel(z);
      [z,iord]=sort(z);
      %     z=[z; 2*z(lz)-z(lz-1)];
      z = [zref; z];
      integrand=[1./H(iord)];
      integrand=[integrand; integrand(lz)];
%        redheight=intrap(integrand,z);
      redheight=cumtrapz(z,integrand);
      redheight=redheight(2:end);
      n1top=n0*exp(-redheight);
      n1sort=zeros(lz,1);
      for iz=1:lz
          n1sort(iord(iz))=n1top(iz);
      end
      ns(inds0,ix2,ix3,1)=n1sort;
    end


    %N+
    ns(:,ix2,ix3,5)=1e-4*ns(:,ix2,ix3,1);

    inds2=inds;
    inds1=setdiff(1:lx1,inds2);


    %MOLECULAR DENSITIES
    nmolc=zeros(lx1,1);
    nmolc(inds1)=(1- rho(inds1)).*ne(inds1);
    if ~isempty(inds2)
      if xg.r(1,1)>xg.r(2,1)
          iref=inds1(1);
      else
        iref=inds1(numel(inds1));
      end
      n0=nmolc(iref);
      ms=30.5*amu;
      H=kb*Tn(inds2,ix2,ix3)./ms./g(inds2,ix2,ix3);
      z = alt(inds2,ix2,ix3);
      lz=numel(z);
      [z,iord]=sort(z);
      z = [z; 2*z(lz)-z(lz-1)];
      integrand=[1./H(iord)];
      integrand=[integrand; integrand(lz)];
%        redheight=intrap(integrand,z);
      redheight=cumtrapz(z,integrand);
      redheight=redheight(2:end);
      nmolctop=n0*exp(-redheight);
      nmolcsort=zeros(lz,1);
      for iz=1:lz
        nmolcsort(iord(iz))=nmolctop(iz);
      end
      nmolc(inds2)=nmolcsort;
    end
    ns(:,ix2,ix3,2) = 1/3 * nmolc;
    ns(:,ix2,ix3,3) = 1/3 * nmolc;
    ns(:,ix2,ix3,4) = 1/3 * nmolc;

    %% PROTONS
    ns(inds2,ix2,ix3,6) = (1 - rho(inds2)) .* ne(inds2);
    z=alt(inds1,ix2,ix3);
    if ~isempty(inds2)
        if xg.r(1,1) > xg.r(2,1)
            iref=inds2(numel(inds2));
        else
            iref=inds2(1);
        end
        n0=ns(iref,ix2,ix3,6);
    else
        [~,iref]=max(alt(:,ix2,ix3));
        n0=1e6;
    end
    ns(inds1,ix2,ix3,6)=chapmana(z,n0,alt(iref,ix2,ix3),mean(Hf));
  end
end
ns(:,:,:,1:6)=max(ns(:,:,:,1:6),mindens);
ns(:,:,:,7)=sum(ns(:,:,:,1:6),4);

vsx1=zeros(lx1,lx2,lx3,7);
Ts=repmat(Tn,[1,1,1,7]);

if closeddip
  % closed dipole grid
  %FIXME:  This code only works for symmetric grids...
  if 2*lx1==xg.lx(1)
    ns=cat(1,ns,ns(lx1:-1:1,:,:,:));
    Ts=cat(1,Ts,Ts(lx1:-1:1,:,:,:));
    vsx1=cat(1,vsx1,vsx1(lx1:-1:1,:,:,:));
  else
    ns=cat(1,ns,ns(lx1,:,:,:),ns(lx1:-1:1,:,:,:));
    Ts=cat(1,Ts,Ts(lx1,:,:,:),Ts(lx1:-1:1,:,:,:));
    vsx1=cat(1,vsx1,vsx1(lx1,:,:,:),vsx1(lx1:-1:1,:,:,:));
  end
end

end
