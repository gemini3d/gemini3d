%GRID RESOLUTIONS
lx1=xg.lx(1);
lx2=xg.lx(2);
lx3=xg.lx(3);

dl1=repmat(xg.dx1f,[1 lx2+4 lx3+4]).*xg.h1;
dl2=repmat(xg.dx2f,[lx1+4 1 lx3+4]).*xg.h2;
dl3=repmat(reshape(xg.dx3f,[1 1 lx3+4]),[lx1+4 lx2+4 1]).*xg.h3;

figure;
imagesc(dl1(:,:,1)/1e3);
colorbar;
title('Resolution along field line (min x3)')
datacursormode;

figure;
imagesc(xg.alt(:,:,1)/1e3);
colorbar;
title('altitude (min x3)')
datacursormode;

figure;
imagesc(dl2(:,:,1)/1e3);
colorbar;
title('Resolution perp. to field line (min x3)')
datacursormode;

figure;
imagesc(squeeze(dl3(1,:,:))/1e3);
colorbar;
title('ZONAL Resolution perp. to field line (min x1)')
datacursormode;

figure;
imagesc(squeeze(dl3(floor(end/2),:,:))/1e3);
colorbar;
title('ZONAL Resolution perp. to field line (mid x1)')
datacursormode;
