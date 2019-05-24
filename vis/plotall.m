function xg = plotall(direc, saveplots, plotfun, xg, visible)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])

narginchk(1,3)
validateattributes(direc, {'char'}, {'vector'}, mfilename, 'path to data', 1)

if nargin<2, saveplots={}; end  %'png', 'eps' or {'png', 'eps'}

if nargin<3
  plotfun=[]; 
else
  validateattributes(plotfun, {'char', 'function_handle'}, {'nonempty'}, mfilename, 'plotting function',3)
end

if nargin<4
  xg=[]; 
else
  validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 4)
end

if nargin<5
  visible = 'on'; 
else
  validateattributes(visible, {'char'}, {'vector'}, mfilename, 'plot visibility: on/off', 5)
end


%% NEED TO READ INPUT FILE TO GET DURATION OF SIMULATION AND START TIME
[ymd0,UTsec0,tdur,dtout]=readconfig([direc,filesep,'inputs/config.ini']);


%% CHECK WHETHER WE NEED TO RELOAD THE GRID (check if one is given because this can take a long time)
if isempty(xg)
  disp('Reloading grid...')
  xg = readgrid([direc,filesep,'inputs',filesep]);
end

plotfun = grid2plotfun(plotfun, xg);

%% TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;
Nt=numel(times);

%% MAIN FIGURE LOOP
% NOTE: keep figure() calls in case plotfcn misses a graphics handle, and
% for Octave...
ymd(1,:) = ymd0;
UTsec(1) = UTsec0;
for i = 2:Nt
  [ymd(i,:), UTsec(i)] = dateinc(dtout, ymd(i-1,:), UTsec(i-1)); %#ok<AGROW>
end

h = plotinit(xg, visible);

if ~isempty(saveplots)
	
  if(isempty(gcp('nocreate')))
    parpool(2);
  end
      
  parfor i = 1:Nt
    plotframe(direc,ymd(i,:),UTsec(i),saveplots,plotfun,xg,h);
  end
else
   
  for i = 1:Nt
    xg = plotframe(direc,ymd(i,:),UTsec(i),saveplots,plotfun,xg,h);
    pause
  end
  
  
end % if saveplots

%% Don't print
if nargout==0, clear('xg'), end
    
end % function

