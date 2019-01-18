function slice3left(hf, x, y, dat, P)

ax = subplot(1,3,1, 'parent', hf, 'nextplot','add','FontSize', P.FS);
%% image
hi = imagesc(ax, x, y, dat);
try
  set(hi, 'alphadata',~isnan(dat));
end
%% line annotation
plot(ax,[P.minxp, P.maxxp],[P.altref, P.altref],'w--','LineWidth',2);
if ~isempty(P.sourcemlat)
  plot(ax,P.sourcemlat,0,'r^','MarkerSize',12,'LineWidth',2);
end
%% axes
slice3axes(ax, P)

xlabel(ax,'eastward dist. (km)');
ylabel(ax,'altitude (km)');

t = datenum(P.ymd(1), P.ymd(2), P.ymd(3), 0, 0, P.utsec);
ttxt = {datestr(t,1), [datestr(t,13),' UT']};
title(ax, ttxt)

%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',18,'Color',[0.66 0.66 0.66],'FontWeight','bold');
%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',16,'Color',[0.5 0.5 0.5],'FontWeight','bold');

end