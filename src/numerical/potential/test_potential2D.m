function test_potential2D(filename)
narginchk(1,1)

addpath([fileparts(mfilename('fullpath')), '/../../../matlab'])

exist_or_skip(filename, 'file')

if isoctave
h = load(filename);

x2 = h.x2;
x3 = h.x3;
Phi = h.Phi;
Phi2 = h.Phi2squeeze;
Phitrue = h.Phitrue;
else
x2 = h5read('/x2');
x3 = h5read('/x3');
Phi = h5read('/Phi');
Phi2 = h5read('/Phi2squeeze');
Phitrue = h5read('/Phitrue');
end


assert_allclose(Phi2(13, 13), 0.00032659, 1e-5,[],'Potential 2d accuracy')

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

% fid=fopen(filename);
% data=fscanf(fid,'%f',1);
% lx2=data(1);
% x2=fscanf(fid,'%f',lx2);
% data=fscanf(fid,'%f',1);
% lx3=data(1);
% x3=fscanf(fid,'%f',lx3);
% Phi=fscanf(fid,'%f',lx2*lx3);
% Phi=reshape(Phi,[lx2,lx3]);
% Phi2=fscanf(fid,'%f',lx2*lx3);
% Phi2=reshape(Phi2,[lx2,lx3]);
% Phitrue=fscanf(fid,'%f',lx2*lx3);
% Phitrue=reshape(Phitrue,[lx2,lx3]);
% fclose(fid);