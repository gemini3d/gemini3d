function [nsi,vs1i,Tsi] = model_resample(xgin,ns,vs1,Ts,xg)
narginchk(5, 5)
validateattributes(xgin, {'struct'}, {'scalar'})
validateattributes(ns, {'numeric'}, {'positive', 'finite'})
validateattributes(vs1, {'numeric'}, {'finite'})
validateattributes(Ts, {'numeric'}, {'positive','finite'})
validateattributes(xg, {'struct'}, {'scalar'})

%% NEW GRID SIZES
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
lsp=size(ns,4);

%% ALLOCATIONS
nsi=zeros(lx1,lx2,lx3,lsp);
vs1i=zeros(lx1,lx2,lx3,lsp);
Tsi=zeros(lx1,lx2,lx3,lsp);

%% INTERPOLATE ONTO NEWER GRID
if lx3 > 1 && lx2 > 1 % 3-D
  disp('interpolating grid for 3-D simulation')
  [X2,X1,X3] = meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2),xgin.x3(3:end-2));
  [X2i,X1i,X3i] = meshgrid(xg.x2(3:end-2),xg.x1(3:end-2),xg.x3(3:end-2));
  for i=1:lsp
    tmpvar = interp3(X2,X1,X3,ns(:,:,:, i),X2i,X1i,X3i);
%      inds=find(isnan(tmpvar));   %this doesn't need to be fixed
%      automatically; the user needs to be aware that they are doing
%      something that is going to likely make the simulation fail or give
%      non-useful results...
%      tmpvar(inds)=1e0;
    nsi(:,:,:, i) = tmpvar;

    tmpvar = interp3(X2,X1,X3,vs1(:,:,:, i),X2i,X1i,X3i);
    vs1i(:,:,:, i) = tmpvar;

    tmpvar = interp3(X2,X1,X3,Ts(:,:,:, i),X2i,X1i,X3i);
    Tsi(:,:,:, i) = tmpvar;
  end
elseif lx3 == 1 % 2-D east-west
  disp('interpolating grid for 2-D simulation in x1, x2')
  [X2,X1]=meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2));
  [X2i,X1i]=meshgrid(xg.x2(3:end-2),xg.x1(3:end-2));
  for i = 1:lsp
    tmpvar=interp2(X2,X1,squeeze(ns(:,:,:, i)),X2i,X1i);
    nsi(:,:,:, i)=reshape(tmpvar,[lx1,lx2,1]);

    tmpvar=interp2(X2,X1,squeeze(vs1(:,:,:, i)),X2i,X1i);
    vs1i(:,:,:, i)=reshape(tmpvar,[lx1,lx2,1]);

    tmpvar=interp2(X2,X1,squeeze(Ts(:,:,:, i)),X2i,X1i);
    Tsi(:,:,:, i)=reshape(tmpvar,[lx1,lx2,1]);
  end
elseif lx2 == 1 % 2-D north-south
  disp('interpolating grid for 2-D simulation in x1, x3')
  % original grid, a priori the first 2 and last 2 values are ghost cells
  % on each axis
  %
  % Detect old non-padded grid and workaround
  if allclose(xgin.x3(1), xg.x3(3), [], 1)
    % old sim, no external ghost cells.
    % Instead of discarding good cells,keep them and say there are
    % new ghost cells outside the grid
    x3in_noghost = linspace(xgin.x3(1), xgin.x3(end), xgin.lx(3));
  else
    % new sim, external ghost cells
    x3in_noghost = xgin.x3(3:end-2);
  end
  x1in_noghost = xgin.x1(3:end-2);
  [X1, X3] = ndgrid(x1in_noghost, x3in_noghost);
  % new grid
  x3_noghost = xg.x3(3:end-2);
  x1_noghost = xg.x1(3:end-2);
  [X1i, X3i] = ndgrid(x1_noghost, x3_noghost);
  % remove degenerate dimension
  ns = permute(ns, [1, 3, 4, 2]);
  vs1 = permute(vs1, [1, 3, 4, 2]);
  Ts = permute(Ts, [1, 3, 4, 2]);
  % for each species
  for i = 1:lsp
    F = griddedInterpolant(X1, X3, ns(:,:, i), 'linear', 'none');
    nsi(:,:,:, i) = F(X1i, X3i);

    F = griddedInterpolant(X1, X3, vs1(:,:, i), 'linear', 'none');
    vs1i(:,:,:, i) = F(X1i, X3i);

    F = griddedInterpolant(X1, X3, Ts(:,:, i), 'linear', 'none');
    Tsi(:,:,:, i) = F(X1i, X3i);
  end
else
  error('Not sure if this is 2-D or 3-D')
end %if

end %function model_resample
