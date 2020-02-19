function plot_grid(path)
%% a high level plot of a 3D grid as a human sanity check
%
% path: path to simgrid.{dat,h5,nc}

xg = readgrid(path);

fig = figure(63);
clf(fig)

ax = subplot(1,3,1, 'parent', fig);
plot(ax, xg.x1/1e3, 'marker','.')
ylabel(ax, 'x1 [km]')
xlabel(ax, 'index (dimensionless)')
title(ax, ['x1 (upward)  N=',int2str(length(xg.x1))])

ax = subplot(1,3,2, 'parent', fig);
plot(ax, xg.x2/1e3, 1:length(xg.x2), 'marker', '.')
xlabel(ax, 'x2 [km]')
ylabel(ax, 'index (dimensionless)')
title(ax, ['x2 (eastward)  N=',int2str(length(xg.x2))])

ax = subplot(1,3,3, 'parent', fig);
plot(ax, xg.x3/1e3, 'marker', '.')
ylabel(ax, 'x3 [km]')
xlabel(ax, 'index (dimensionless)')
title(ax, ['x3 (northward)  N=',int2str(length(xg.x3))])
end
