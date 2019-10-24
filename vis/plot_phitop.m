function plot_phitop(x, y, Phitop, h, P)
narginchk(5, 5)
validateattr(x, {'numeric'}, {'vector'}, mfilename, 'x distance', 1)
validateattr(y, {'numeric'}, {'vector'}, mfilename, 'y distance', 2)
validateattr(Phitop, {'numeric'}, {'2d'}, mfilename, 'precipitation', 3)

ax = get_axes(h);
hi = imagesc(x, y, Phitop, 'parent', ax);
try %#ok<TRYNC> octave < 5
  set(hi, 'alphadata', ~isnan(dat));
end


if ndims(Phitop) == 1

else
  imagesc(Phitop, 'parent', ax)
  colorbar('peer', ax)
  ylabel(ax, 'northward dist. (km)');
  xlabel(ax, 'eastward dist. (km)');
  slice3axes(ax, P)
end