function test_potential2D(filename)
narginchk(1,1)

cwd = fileparts(mfilename('fullpath'));
addpath([cwd,filesep,'..',filesep,'..', filesep, 'script_utils'])
addpath([cwd,filesep,'..',filesep,'..', filesep, 'tests'])

fid=fopen(filename);
data=fscanf(fid,'%f',1);
lx2=data(1);
x2=fscanf(fid,'%f',lx2);
data=fscanf(fid,'%f',1);
lx3=data(1);
x3=fscanf(fid,'%f',lx3);
Phi=fscanf(fid,'%f',lx2*lx3);
Phi=reshape(Phi,[lx2,lx3]);
Phi2=fscanf(fid,'%f',lx2*lx3);
Phi2=reshape(Phi2,[lx2,lx3]);
Phitrue=fscanf(fid,'%f',lx2*lx3);
Phitrue=reshape(Phitrue,[lx2,lx3]);
fclose(fid);

assert_allclose(Phi2(13, 13), 0.000327, 1e-5,[],'Potential 2d accuracy')

if ~isinteractive
  return
end

%% Plot data
figure(1);

subplot(1,3,1)
imagesc(x2,x3,Phi);
colorbar;
axis xy;
xlabel('distance (m)')
ylabel('distance (m)')
title('2D potential (polarization)')

subplot(1,3,2)
imagesc(x2,x3,Phi2);
colorbar;
axis xy;
xlabel('distance (m)')
ylabel('distance (m)')
title('2D potential (static)')

subplot(1,3,3)
imagesc(x2,x3,Phitrue);
colorbar;
axis xy;
xlabel('distance (m)')
ylabel('distance (m)')
title('2D potential (analytical)')

end