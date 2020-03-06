function plot2D_cart(ymd,UTsec,xg,parm,parmlbl,caxlims,sourceloc,ha, cmap)

narginchk(4,9)
validateattributes(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 1)
validateattributes(UTsec, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 2)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 3)
validateattributes(parm, {'numeric'}, {'real'}, mfilename, 'parameter to plot',4)
if nargin<5, parmlbl=''; end
validateattributes(parmlbl, {'char'}, {'vector'}, mfilename, 'parameter label', 5)

if nargin<6
  caxlims=[];
else
  validateattributes(caxlims, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'plot intensity (min, max)', 6)
end
if nargin<7  || isempty(sourceloc) % leave || for validate
  sourceloc = [];
else
  validateattributes(sourceloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'source magnetic coordinates', 7)
end

if nargin<8 || isempty(ha)
  ha = get_axes();
else
  ha = get_axes(ha);
end

if nargin<9 || isempty(cmap)
  cmap = parula(256);
end


%% PLOT parameters
%set(h,'PaperPosition',[0 0 11 4.5]);


%SOURCE LOCATION
if ~isempty(sourceloc)
  sourcemlat=sourceloc(1);
  %sourcemlon=sourceloc(2);
else
  sourcemlat=[];
  %sourcemlon=[];
end


%SIZE OF SIMULATION
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
inds1=3:lx1+2;
inds2=3:lx2+2;
inds3=3:lx3+2;
Re=6370e3;


%JUST PICK AN X3 LOCATION FOR THE MERIDIONAL SLICE PLOT, AND AN ALTITUDE FOR THE LAT./LON. SLICE
%ix3=floor(lx3/2);
altref=300;


%SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
meantheta=mean(xg.theta(:));
%meanphi=mean(xg.phi(:));
y=-1*(xg.theta-meantheta);   %this is a mag colat. coordinate and is only used for defining grid in linspaces below, runs backward from north distance, hence the negative sign
%x=(xg.phi-meanphi);       %mag. lon coordinate, pos. eastward
x=xg.x2(inds2)/Re/sin(meantheta);
z=xg.alt;
lxp=500;
lyp=500;
lzp=500;
minx=min(x(:));
maxx=max(x(:));
miny=min(y(:));
maxy=max(y(:));
minz=min(z(:));
maxz=max(z(:));
xp=linspace(minx,maxx,lxp);     %eastward distance (rads.)
yp=linspace(miny,maxy,lyp);     %should be interpreted as northward distance (in rads.).  Irrespective of ordering of xg.theta, this will be monotonic increasing!!!
zp=linspace(minz,maxz,lzp)';     %altitude (meters)


assert(ismatrix(parm) && ~isscalar(parm), ['need 2-D or 1-D ',parmlbl,' -- is squeeze() needed?  size(parm):', num2str(size(parm))])

%INTERPOLATE ONTO PLOTTING GRID
if xg.lx(3)==1     %alt./lon. slice
  if isvector(parm)
    parmp = interp1(xg.x2(inds2), parm, xp*Re*sin(meantheta));
  elseif ismatrix(parm)
    [X,Z]=meshgrid(xp,zp);    %meridional meshgrid, this defines the grid for plotting
    x1plot=Z(:);   % upward distance (meters)
    x2plot=X(:)*Re*sin(meantheta);     % eastward distance (meters)

    parmtmp=parm(:,:);
    parmp=interp2(xg.x2(inds2),xg.x1(inds1),parmtmp,x2plot,x1plot);
    parmp=reshape(parmp,[lzp,lxp]);    %slice expects the first dim. to be "y" ("z" in the 2D case)
  end
elseif xg.lx(2)==1    %alt./lat. slice
  if isvector(parm)
    parmp = interp1(xg.x3(inds3), parm, yp*Re);
  elseif ismatrix(parm)
    [Y,Z]=meshgrid(yp,zp);

    x1plot = Z(:); % upward distance (meters)
    x3plot = Y(:)*Re; % northward distance (meters)

    %ix2=floor(lx2/2);
    parmtmp=parm(:,:);     %so north dist, east dist., alt.
    parmp=interp2(xg.x3(inds3),xg.x1(inds1),parmtmp,x3plot,x1plot);
    parmp=reshape(parmp,[lzp,lyp]);    %slice expects the first dim. to be "y"
  end
end


%CONVERT ANGULAR COORDINATES TO MLAT,MLON
if (xg.lx(2)==1)
  yp=yp*Re;
  [yp,inds]=sort(yp);
  %parmp2=parmp2(inds,:,:);
  parmp=parmp(:,inds);
elseif (xg.lx(3)==1)
  %xp=(xp+meanphi)*180/pi;
  xp=xp*Re*sin(meantheta);    %eastward ground distance
  [xp,inds]=sort(xp);
  parmp=parmp(:,inds,:);
  %parmp2=parmp2(:,inds,:);
end

%% BOUNDS FOR THE PLOTTING
minxp=min(xp(:));
maxxp=max(xp(:));
%minyp=min(yp(:));
%maxyp=max(yp(:));
%minzp=min(zp(:));
%maxzp=max(zp(:));

%% MAKE THE PLOT!
if isvector(parm)
  if xg.lx(3) == 1
    plot1d2(xp, parmp, ha, parmlbl)
  elseif xg.lx(2) == 1
    plot1d3(yp, parmp, ha, parmlbl)
  end
else
  if xg.lx(3)==1
    plot12(xp, zp, parmp, ha, sourcemlat, minxp, maxxp, altref, cmap, caxlims, parmlbl)
  elseif xg.lx(2)==1
    plot13(yp, zp, parmp, ha, sourcemlat, cmap, caxlims, parmlbl)
  end
end
ttxt = time2str(ymd, UTsec);
title(ha, ttxt)
%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',18,'Color',[0.66 0.66 0.66],'FontWeight','bold');
%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',16,'Color',[0.5 0.5 0.5],'FontWeight','bold');

end

function plot1d2(xp, parm, ha, parmlbl)
narginchk(4,4)
plot(ha, xp/1e3, parm)
xlabel(ha, 'eastward dist. (km)')
ylabel(ha, parmlbl)
end % function


function plot1d3(yp, parm, ha, parmlbl)
narginchk(4,4)
plot(ha, yp/1e3, parm)
xlabel(ha, 'northward dist. (km)')
ylabel(ha, parmlbl)
end % function


function plot12(xp, zp, parmp, ha, sourcemlat, minxp, maxxp, altref, cmap, caxlims, parmlbl)
narginchk(11, 11)
hi = imagesc(xp/1e3, zp/1e3,parmp, 'parent', ha);
hold(ha, 'on')
plot(ha, [minxp/1e3,maxxp/1e3],[altref, altref],'w--','LineWidth',2)

if ~isempty(sourcemlat)
  plot(ha, sourcemlat,0,'r^','MarkerSize',12,'LineWidth',2)
end
hold(ha, 'off')

try %#ok<*TRYNC> % octave < 5
  set(hi, 'alphadata', ~isnan(parmp))
end

make_colorbar(ha, cmap, caxlims, parmlbl)

xlabel(ha, 'eastward dist. (km)')
ylabel(ha, 'altitude (km)')

end

function plot13(yp, zp, parmp, ha, sourcemlat, cmap, caxlims, parmlbl)
narginchk(8,8)
hi = imagesc(yp/1e3, zp/1e3, parmp, 'parent', ha);
hold(ha, 'on')
%plot([minyp,maxyp],[altref,altref],'w--','LineWidth',2);
if (~isempty(sourcemlat))
  plot(ha, sourcemlat,0,'r^','MarkerSize',12,'LineWidth',2);
end
hold(ha, 'off')

try % % octave < 5
  set(hi, 'alphadata', ~isnan(parmp));
end

make_colorbar(ha, cmap, caxlims, parmlbl)

xlabel(ha, 'northward dist. (km)')
ylabel(ha, 'altitude (km)')

end % function


function make_colorbar(ha, cmap, lims, parmlbl)

tight_axis(ha)
colormap(ha, cmap)
if ~isempty(lims) && all(~isnan(lims)) && lims(1) < lims(2)
  % keep semicolon here as caxis can be buggy and output to console
  caxis(ha, lims);
end
c=colorbar('peer', ha);
xlabel(c, parmlbl)

end % function