function slice3axes(ax, P)

tight_axis(ax)
colormap(ax, P.cmap)
caxis(ax, P.caxlims);

c = colorbar(ax);
xlabel(c, P.parmlbl);

end