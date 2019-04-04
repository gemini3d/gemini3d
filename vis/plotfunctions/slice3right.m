function slice3right(hf, x, y, dat, P)

ax = subplot(1,3,3, 'parent', hf, 'nextplot','add','FontSize',P.FS);
%% image
hi = imagesc(ax, x, y, dat);
try  % octave < 5
  set(hi, 'alphadata', ~isnan(dat));
end
%% line
% plot(ha, [minyp,maxyp],[altref,altref],'w--','LineWidth',2);
if ~isempty(P.sourcemlat)
  plot(ax, P.sourcemlat,0,'r^','MarkerSize',12,'LineWidth',2);
end
%% axes
slice3axes(ax, P)

xlabel(ax, 'northward dist. (km)');
ylabel(ax, 'altitude (km)');
end