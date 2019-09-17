function [x1,x2,x3,f]=testinterp3(filename, realbits)

narginchk(1,2)
cwd = fileparts(mfilename('fullpath'));
addpath([cwd,filesep,'..',filesep,'..',filesep,'script_utils'])

if nargin == 1
  realbits = 64;
end

validateattr(realbits, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename,'real bits',2)

exist_or_skip(filename, 'file')

switch realbits
  case 64, freal = 'float64';
  case 32, freal = 'float32';
  otherwise, error(['unknown precision', num2str(realbits)])
end

fid=fopen(filename,'r');

%% LOAD DATA
lx1=fread(fid,1,'integer*4');
lx2=fread(fid,1,'integer*4');
lx3=fread(fid,1,'integer*4');
x1=fread(fid,lx1, freal);
x2=fread(fid,lx2, freal);
x3=fread(fid,lx3, freal);
f=fread(fid,lx1*lx2*lx3, freal);
f=reshape(f,[lx1,lx2,lx3]);

fclose(fid);

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
