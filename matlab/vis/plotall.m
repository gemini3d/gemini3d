function xg = plotall(direc, saveplot_fmt, plotfun, xg, visible)

narginchk(1,5)

validateattributes(direc, {'char'}, {'vector'}, mfilename, 'path to data', 1)
direc = absolute_path(direc);
assert(is_folder(direc), [direc, ' is not a directory'])

if nargin<2, saveplot_fmt={}; end  %e.g. {'png'} or {'png', 'eps'}

if nargin<3, plotfun=[]; end
if ~isempty(plotfun)
  validateattributes(plotfun, {'char', 'function_handle'}, {'nonempty'}, mfilename, 'plotting function',3)
end

if nargin<4, xg=[]; end
if ~isempty(xg)
  validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 4)
end

if nargin<5
  if isempty(saveplot_fmt)
    visible = 'on';
  else
    visible = 'off';
  end
end
validateattributes(visible, {'char'}, {'vector'}, mfilename, 'plot visibility: on/off', 5)

lxs = simsize(direc);
disp(['sim grid dimensions: ',num2str(lxs)])


%% NEED TO READ INPUT FILE TO GET DURATION OF SIMULATION AND START TIME
params = read_config([direc, '/inputs']);

%% CHECK WHETHER WE NEED TO RELOAD THE GRID (check if one is given because this can take a long time)
if isempty(xg)
  disp('Reloading grid...')
  xg = readgrid([direc, '/inputs']);
end

plotfun = grid2plotfun(plotfun, xg);

%% TIMES OF INTEREST
times=params.UTsec0:params.dtout:params.UTsec0+params.tdur;
Nt=numel(times);

%% MAIN FIGURE LOOP
% NOTE: keep figure() calls in case plotfcn misses a graphics handle, and
% for Octave...
ymd(1,:) = params.ymd;
UTsec(1) = params.UTsec0;
for i = 2:Nt
  [ymd(i,:), UTsec(i)] = dateinc(params.dtout, ymd(i-1,:), UTsec(i-1)); %#ok<AGROW>
end

h = plotinit(xg, visible);

if ~isempty(saveplot_fmt)
  % plot and save as fast as possible.
  for i = 1:Nt
    plotframe(direc, ymd(i,:), UTsec(i), saveplot_fmt, plotfun, xg, h)
  end
elseif isoctave || isinteractive
  % displaying interactively, not saving
  % Note: CLI Octave can plot also
  for i = 1:Nt
    plotframe(direc, ymd(i,:), UTsec(i), saveplot_fmt, plotfun, xg, h)

    drawnow % need this here to ensure plots update (race condition)
    disp(''), disp('** press any key to plot next time step, or Ctrl C to stop**')
    pause
  end
else
  error('No Matlab / Octave desktop so cannot plot. Was also not told to save')
end % if saveplots

if is_folder([direc, '/aurmaps']) % glow sim
  plotglow(direc, saveplot_fmt, visible)
end

%% Don't print
if nargout==0, clear('xg'), end

end % function
