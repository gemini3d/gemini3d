%LOAD AND PLOT NUMERICAL SOLUTION
function test_diffusion1D(fn)
narginchk(1,1)
cwd = fileparts(mfilename('fullpath'));
addpath([cwd,filesep,'../../tests'])
addpath([cwd,filesep,'../../script_utils'])

exist_or_skip(fn)

fid=fopen(fn);
data=fscanf(fid,'%f',2);
lt=data(1);
lx1=data(2);
x1=fscanf(fid,'%f',lx1+4)';

%M=fscanf(fid,'%f',lx1*5)';
%M=reshape(M,[lx1 5])';
%b=fscanf(fid,'%f',lx1)';

Ts=zeros(lx1,lt);
t=zeros(lt,1);
for it=1:lt
  t(it)=fscanf(fid,'%f',1);
  Ts(:,it)=fscanf(fid,'%f',lx1)';
end % for

% reltol = 1e-5 for real32
assert_allclose(Ts(13,end), 0.2757552094055,1e-5,[],'1-D diffusion accuracy')

if ~isinteractive
  return
end
%% plots
figure
imagesc(t,x1(3:end-2),Ts)
colorbar
xlabel('time [sec]')
ylabel('displacement [m]')
title('1-D diffusion (vs. time)')

end % function
