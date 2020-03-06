function plotglow(direc, saveplot_fmt, visible)
% plots Gemini-Glow auroral emissions
narginchk(1,3)
if nargin<2, saveplot_fmt={'png'}; end  %e.g. {'png'} or {'png', 'eps'}

assert(is_folder(direc), [direc, ' is not a directory.'])
config_dir = [direc, '/inputs'];
aurora_dir = [direc, '/aurmaps'];

%array of volume emission rates at each altitude; cm-3 s-1:
wavelengths = {'3371', '4278', '5200', '5577', '6300', '7320', '10400', ...
  '3466', '7774', '8446', '3726', 'LBH', '1356', '1493', '1304'};

%READ IN SIMULATION INFO
params = read_config(config_dir);
%glow_params = read_nml_group([config_dir, '/config.nml'], 'glow');

%READ IN THE GRID
xg = readgrid(config_dir);

%% GET THE SYSTEM SIZE
lwave=length(wavelengths);
lx2=xg.lx(2);
lx3=xg.lx(3);
x2=xg.x2(3:end-2);
x3=xg.x3(3:end-2);

%% get file list
for ext = {'.h5', '.dat'}
  file_list = dir([aurora_dir, '/*', ext{1}]);
  if ~isempty(file_list), break, end
end
assert(~isempty(file_list), ['No auroral data found in ', aurora_dir])

%% make plots
ymd = params.ymd;
UTsec = params.UTsec0;
time_str = time2str(ymd, UTsec);
hf = [];
for i = 1:length(file_list)
  filename = [aurora_dir, '/', file_list(i).name];
  bFrame = loadglow_aurmap(filename, lx2, lx3, lwave);

  if lx3 > 1  % 3D sim
    hf = plot_emission_line(x2, x3, bFrame, time_str, wavelengths, hf, visible);
  else  % 2D sim
    hf = plot_emissions(x2, wavelengths, squeeze(bFrame), time_str, hf, visible);
  end

  save_glowframe(params.flagoutput, filename, saveplot_fmt, hf)

  % we use dtout instead of dtglow because we're only plotting times the
  % full simulation outputs too.
  [ymd, UTsec] = dateinc(params.dtout, ymd, UTsec);
  time_str = time2str(ymd, UTsec);
end

end % function


function hf = plot_emissions(x2, wavelengths, bFrame, time_str, hf, visible)
narginchk(4,6)
validateattributes(x2, {'numeric'}, {'vector'}, mfilename)
validateattributes(bFrame, {'numeric'}, {'ndims', 2}, mfilename)
validateattributes(wavelengths, {'cell'}, {'vector'}, mfilename)

if nargin < 5 || isempty(hf)
  hf = make_glowfig(visible);
else
  clf(hf)
end

ax = axes('parent', hf);
imagesc(1:length(wavelengths), x2/1e3,squeeze(bFrame))    % singleton dimension since 2D simulation
set(ax, 'xtick', 1:length(wavelengths), 'xticklabel', wavelengths)
ylabel(ax, 'Eastward Distance (km)')
xlabel(ax, 'emission wavelength (\AA)', 'interpreter', 'latex')
title(ax, time_str)
hc = colorbar('peer', ax);
ylabel(hc, 'Intensity (R)')
end


function hf = plot_emission_line(x2, x3, bFrame, time_str, wavelengths, hf, visible)
narginchk(5,7)
validateattributes(x2, {'numeric'}, {'vector'}, mfilename)
validateattributes(x3, {'numeric'}, {'vector'}, mfilename)
validateattributes(bFrame, {'numeric'}, {'ndims', 3}, mfilename)
validateattributes(time_str, {'cell'}, {'vector'}, mfilename)
validateattributes(wavelengths, {'cell'}, {'vector'}, mfilename)

if nargin < 5 || isempty(hf)
  hf = make_glowfig(visible);
else
  clf(hf)
end

% arbitrary pick of which emission lines to plot lat/lon slices
inds = [2, 4, 5, 9];

for i=1:length(inds)
  ax = subplot(2,2,i,'parent', hf);
  imagesc(x2/1e3, x3/1e3, squeeze(bFrame(:,:,inds(i)))', 'parent', ax);
  axis(ax, 'xy')
  axis(ax, 'tight')
  %caxis(caxlims);
  cb = colorbar('peer', ax);
  %set(cb,'yticklabel',sprintf('10^{%g}|', get(cb,'ytick')))
  ylabel(cb,'Intensity (R)')
  title(ax, [wavelengths{inds(i)},'\AA  intensity: ', time_str{1}, ' ', time_str{2}], 'interpreter', 'latex')
end

ax=subplot(2,2,3);
xlabel(ax, 'Eastward Distance (km)')
ylabel(ax, 'Northward Distance (km)')

end % function


function hf = make_glowfig(visible)
narginchk(0,1)
if nargin < 1, visible = 'on'; end
validateattributes(visible, {'char'}, {'vector'})

hf = figure('toolbar', 'none');
pos = get(hf, 'position');
set(hf, 'name', 'aurora', 'unit', 'pixels', 'position', [pos(1), pos(2), 800, 500], 'visible', visible)

end
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libxvid -r 30 -q:v 0 isinglass_geminiglow_4278.avi
%ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
