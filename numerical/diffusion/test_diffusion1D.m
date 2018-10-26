%LOAD AND PLOT NUMERICAL SOLUTION
function test_diffusion1D(fn)

  cwd = fileparts(mfilename('fullpath'));
  addpath([cwd,filesep,'../../tests'])

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
  %  plot(x1(3:end-2),Ts);
  %  pause(0.1)
  end % for

  figure
  imagesc(t,x1(3:end-2),Ts)
  colorbar
  xlabel('time [sec]')
  ylabel('displacement [m]')
  title('1-D diffusion (vs. time)')
  
  assert_allclose(Ts(13,end),0.2757552094055,[],[],'1-D diffusion accuracy')
end % function

