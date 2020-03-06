function [x1,x2,x3,f] = testinterp3(filename)

narginchk(1,1)

addpath([fileparts(mfilename('fullpath')), '/../../../matlab'])

exist_or_skip(filename, 'file')

if isoctave
  h = load(filename);
  lx1 = h.lx1;
  lx2 = h.lx2;
  lx3 = h.lx3;
  x1 = h.x1;
  x2 = h.x2;
  x3 = h.x3;
  f = h.f;
else
  lx1 = h5read(filename, '/lx1');
  lx2 = h5read(filename, '/lx2');
  lx3 = h5read(filename, '/lx3');
  x1 = h5read(filename, '/x1');
  x2 = h5read(filename, '/x2');
  x3 = h5read(filename, '/x3');
  f = h5read(filename, '/f');
end

assert(lx1==256, int2str(size(lx1)))
assert(lx2==256, int2str(size(lx2)))
assert(lx3==256, int2str(size(lx3)))
assert(all(size(f) == [lx1,lx2,lx3]), 'array size mismatch')

if ~isinteractive
  if ~nargout, clear, end
  return
end

%% PLOT
figure

subplot(1,3,1)
imagesc(x2,x1,f(:,:,end/2));
axis xy;
xlabel('x_2')
ylabel('x_1')
c=colorbar;
ylabel(c,'f')
title('3D interp x_1-x_2')

subplot(1,3,2)
imagesc(x3,x1,squeeze(f(:,end/2-10,:)));
axis xy;
xlabel('x_3')
ylabel('x_1')
c=colorbar;
ylabel(c,'f')
title('3D interp x_1-x_3')

subplot(1,3,3)
imagesc(x2,x3,squeeze(f(end/2-10,:,:)));
axis xy;
xlabel('x_2')
ylabel('x_3')
c=colorbar;
ylabel(c,'f')
title('3D interp x_2-x_3')
end % function

% fid=fopen(filename,'r');

% %% LOAD DATA
% lx1=fread(fid,1,'integer*4');
% lx2=fread(fid,1,'integer*4');
% lx3=fread(fid,1,'integer*4');
% x1=fread(fid,lx1, freal);
% x2=fread(fid,lx2, freal);
% x3=fread(fid,lx3, freal);
% f=fread(fid,lx1*lx2*lx3, freal);
% f=reshape(f,[lx1,lx2,lx3]);

% fclose(fid);