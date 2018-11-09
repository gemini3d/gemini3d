function xg = plotall(direc,saveplots,plotfun,xg)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])

narginchk(1,3)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'path to data', 1)

if nargin<2, saveplots={}; end  %'png', 'eps' or {'png', 'eps'}

if nargin<3, plotfun=[]; end    %need to validate input???

if nargin<4
  xg=[]; 
else
  validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 4)
end


%%NEED TO READ INPUT FILE TO GET DURATION OF SIMULATION AND START TIME
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,filesep,'inputs/config.ini']);


%%CHECK WHETHER WE NEED TO RELOAD THE GRID (check if one is given because this can take a long time)
if (isempty(xg))
  disp('Reloading grid...  Provide one as input if you do not want this to happen...\n')
  xg = readgrid([direc,filesep,'inputs',filesep]);
end


%%DEFINE THE PLOTTING FUNCTION BASED ON THE TYPE OF GRID USED
if (isempty(plotfun))
  minh1=min(xg.h1(:));
  maxh1=max(xg.h1(:));
  if (abs(minh1-1)>1e-4 | abs(maxh1-1)>1e-4)    %curvilinear grid
    if (xg.lx(2)>1 & xg.lx(3)>1)
      plotfun=@plot3D_curv_frames_long;
    else
      plotfun=@plot2D_curv;
    end
  else     %cartesian grid
    if (xg.lx(2)>1 & xg.lx(3)>1)
      plotfun=@plot3D_cart_frames_long_ENU;
    else
      plotfun=@plot2D_cart;
    end 
  end
end


%% TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;
lt=numel(times);


%%INITIALIZE FIGURE SETS
Csp = ceil(sqrt(lt));
Rsp = ceil(lt/Csp);

h.f1=figure('name','V1', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f2=figure('name','Ti', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f3=figure('name','Te', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f4=figure('name','J1', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f5=figure('name','V2', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f6=figure('name','V3', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f7=figure('name','J2', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f8=figure('name','J3', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
if xg.lx(2)>1 && xg.lx(3)>1 % a 3-D simulation
  h.f9=figure('name','phiTop', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
else
  h.f9 = [];
end
h.f10=figure('name','Ne', 'units', 'normalized', 'position', [.1, .1, .5, .5]);

lotsplots = ~isempty(h.f9) || lt > 16;


%%MAIN FIGURE LOOP
% NOTE: keep figure() calls in case plotfcn misses a graphics handle, and
% for Octave...
ymd=ymd0;
UTsec=UTsec0;
for it=1:lt
    xg=plotframe(direc,ymd,UTsec,saveplots,plotfun,xg,h);
    [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end % for


%%UNLOAD GRID DATA IF USER DOESN'T WANT IT 
if nargout==0, clear('xg'), end
    
end % function

