function plot_phitop(x, y, Phitop, h, P)
% phi_top is Topside Potential

narginchk(5, 5)
validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x distance', 1)
validateattributes(y, {'numeric'}, {'vector'}, mfilename, 'y distance', 2)
validateattributes(Phitop, {'numeric'}, {'2d'}, mfilename, 'Potential', 3)
validateattributes(P, {'struct'}, {'scalar'}, mfilename, 'parameters', 5)

ax = get_axes(h);
hi = imagesc(x, y, Phitop, 'parent', ax);
try %#ok<TRYNC> octave < 5
  set(hi, 'alphadata', ~isnan(dat));
end

axes_tidy(ax, P)

ylabel(ax, 'northward dist. (km)');
xlabel(ax, 'eastward dist. (km)');

title(ax, time2str(P.ymd, P.utsec))


end