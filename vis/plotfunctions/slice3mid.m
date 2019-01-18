function slice3mid(hf, x, y, dat, P)

ax = subplot(1,3,2, 'parent', hf, 'nextplot','add','FontSize', P.FS);
%% image
hi = imagesc(ax, x, y, dat);
try
  set(hi, 'alphadata',~isnan(dat));
end
%% line
if ~isempty(P.sourcemlat)
  plot(ax,[P.minxp, P.maxxp],[P.sourcemlon, P.sourcemlon],'w--','LineWidth',2);
  plot(ax, P.sourcemlat, P.sourcemlon,'r^','MarkerSize',12,'LineWidth',2);
end

%% axes
slice3axes(ax, P)

ylabel(ax,'northward dist. (km)');
xlabel(ax,'eastward dist. (km)');

end