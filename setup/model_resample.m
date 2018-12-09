function [nsi,vs1i,Tsi]=model_resample(xgin,ns,vs1,Ts,xg)

  %%NEW GRID SIZES
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  lsp=size(ns,4);


  %$ALLOCATIONS
  nsi=zeros(lx1,lx2,lx3,lsp);
  vs1i=zeros(lx1,lx2,lx3,lsp);
  Tsi=zeros(lx1,lx2,lx3,lsp);


  %%INTERPOLATE ONTO NEWER GRID
  if (lx3~=1)
    fprintf('Starting interp3''s...\n');
    [X2,X1,X3]=meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2),xgin.x3(3:end-2));
    [X2i,X1i,X3i]=meshgrid(xg.x2(3:end-2),xg.x1(3:end-2),xg.x3(3:end-2));
    for isp=1:lsp
      tmpvar=interp3(X2,X1,X3,ns(:,:,:,isp),X2i,X1i,X3i);
%      inds=find(isnan(tmpvar));   %this doesn't need to be fixed
%      automatically; the user needs to be aware that they are doing
%      something that is going to likely make the simulation fail or give
%      non-useful results...
%      tmpvar(inds)=1e0;
      nsi(:,:,:,isp)=tmpvar;
      tmpvar=interp3(X2,X1,X3,vs1(:,:,:,isp),X2i,X1i,X3i);
%      tmpvar(inds)=0e0;
      vs1i(:,:,:,isp)=tmpvar;
      tmpvar=interp3(X2,X1,X3,Ts(:,:,:,isp),X2i,X1i,X3i);
%      tmpvar(inds)=100e0;
      Tsi(:,:,:,isp)=tmpvar;
    end
  else
    fprintf('Starting interp2''s...\n');
    [X2,X1]=meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2));
    [X2i,X1i]=meshgrid(xg.x2(3:end-2),xg.x1(3:end-2));
    for isp=1:lsp
      tmpvar=interp2(X2,X1,squeeze(ns(:,:,:,isp)),X2i,X1i);
%      inds=find(isnan(tmpvar));
%      tmpvar(inds)=1e0;
      nsi(:,:,:,isp)=reshape(tmpvar,[lx1,lx2,1]);
      tmpvar=interp2(X2,X1,squeeze(vs1(:,:,:,isp)),X2i,X1i);
%      tmpvar(inds)=0e0;
      vs1i(:,:,:,isp)=reshape(tmpvar,[lx1,lx2,1]);
      tmpvar=interp2(X2,X1,squeeze(Ts(:,:,:,isp)),X2i,X1i);
%      tmpvar(inds)=100e0;
      Tsi(:,:,:,isp)=reshape(tmpvar,[lx1,lx2,1]);
    end
  end %if

end %function model_resample
