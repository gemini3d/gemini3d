function axes_tidy(ax, P)
narginchk(2,2)

validateattributes(P, {'struct'}, {'scalar'}, mfilename)

tight_axis(ax)
colormap(ax, P.cmap)
caxis(ax, P.caxlims);

c = colorbar(ax);
xlabel(c, P.parmlbl);

end