%LOAD AND PLOT NUMERICAL SOLUTION to diffusion problem
%return the test data to the user in case they want to
%look over the differences and numerical error.
function [x1,TsEuler,TsBDF2,Tstrue]=test_diffusion1D(fn)

narginchk(1,1)
cwd = fileparts(mfilename('fullpath'));
for p = {'tests', 'script_utils', 'vis'}
  addpath([cwd,filesep,'..',filesep,'..', filesep, p{:}])
end

exist_or_skip(fn, 'file')

fid=fopen(fn);
data=fscanf(fid,'%f',2);
lt=data(1);
lx1=data(2);
x1=fscanf(fid,'%f',lx1+4)';

t=zeros(lt,1);
for it=1:lt
  t(it)=fscanf(fid,'%f',1);
  TsEuler(:,it)=fscanf(fid,'%f',lx1)';
  TsBDF2(:,it)=fscanf(fid,'%f',lx1)';
  Tstrue(:,it)=fscanf(fid,'%f',lx1)';
end % for

% reltol = 1e-5 for real32
% this is just a random point we're comparing.
assert_allclose(TsEuler(13, 13), 0.770938954253086, 1e-5,[],'1-D Euler diffusion accuracy')
assert_allclose(TsBDF2(13, 13),  0.763236513549944, 1e-5,[],'1-D BDF2 diffusion accuracy')
assert_allclose(Tstrue(13, 13),  0.763014494788105, 1e-5,[],'1-D true diffusion accuracy')

if ~isinteractive
  if ~nargout, clear, end
  return
end
%% plots

figure
ax1=subplot(1,3,1);
imagesc(t,x1(3:end-2),TsEuler, 'parent', ax1)
ylabel(ax1,'distance (m)')
title(ax1,'1D diffusion (backward Euler)')

ax2=subplot(1,3,2);
imagesc(t,x1(3:end-2),TsBDF2, 'parent', ax2)
title(ax2,'1D diffusion (TRBDF2)')

ax3=subplot(1,3,3);
imagesc(t,x1(3:end-2),Tstrue, 'parent', ax3)
title(ax3,'1D diffusion (analytic)')

for ax = [ax1,ax2,ax3]
  colormap(ax, 'bwr')
  colorbar('peer',ax)
  xlabel(ax,'time (sec)')
end

if ~nargout, clear, end
end % function
