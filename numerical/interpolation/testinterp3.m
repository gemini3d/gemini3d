function [x1,x2,x3,f]=testinterp3(filename)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd,filesep,'..',filesep,'..',filesep,'script_utils'])

validateattr(filename, {'char'}, {'vector'}, mfilename,'interp file to compare',1)

%% if path not exist, exit code 77 as GNU standard skip test indicator
if exist(filename, 'file') ~= 2, fprintf(2,[filename,' not found\n']), exit(77), end

fid=fopen(filename,'r');

%% LOAD DATA
lx1=fread(fid,1,'integer*4');
lx2=fread(fid,1,'integer*4');
lx3=fread(fid,1,'integer*4');
x1=fread(fid,lx1,'real*8');
x2=fread(fid,lx2,'real*8');
x3=fread(fid,lx3,'real*8');
f=fread(fid,lx1*lx2*lx3,'real*8');
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
