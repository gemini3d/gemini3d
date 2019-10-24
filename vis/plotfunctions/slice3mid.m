function slice3mid(hf, x, y, dat, P)
narginchk(5,5)

ax = subplot(1,3,2, 'parent', hf, 'nextplot','add','FontSize', P.FS);
%% image
hi = imagesc(x, y, dat, 'parent', ax);
try %#ok<TRYNC> octave < 5
  set(hi, 'alphadata', ~isnan(dat));
end
%% line
if ~isempty(P.sourcemlat)
  plot(ax,[P.minxp, P.maxxp],[P.sourcemlon, P.sourcemlon],'w--','LineWidth',2);
  plot(ax, P.sourcemlat, P.sourcemlon,'r^','MarkerSize',12,'LineWidth',2);
end

%% axes
slice3axes(ax, P)

ylabel(ax, 'northward dist. (km)');
xlabel(ax, 'eastward dist. (km)');

end
