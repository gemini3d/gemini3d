function dat = loadframe3Dcurv_raw(filename)

narginchk(1,1)
%% SIMULATION SIZE
lsp=7;
lxs = simsize(filename);
%% SIMULATION RESULTS
assert(is_file(filename), [filename,' is not a file.'])
dat.filename = filename;

fid=fopen(filename,'r');
simdt(fid);
%% load densities
dat.ns = read4D(fid, lsp, lxs);
%% load Vparallel
dat.vs1 = read4D(fid, lsp, lxs);
dat.v1 = sum(dat.ns(:,:,:,1:6) .* dat.vs1(:,:,:,1:6),4) ./ dat.ns(:,:,:,lsp);
%% load temperatures
dat.Ts = read4D(fid, lsp, lxs);
%% load current densities
dat.J1 = read3D(fid, lxs);
dat.J2 = read3D(fid, lxs);
dat.J3 = read3D(fid, lxs);
%% load Vperp
dat.v2 = read3D(fid, lxs);
dat.v3 = read3D(fid, lxs);
%% load topside potential
dat.Phitop = read2D(fid, lxs);

fclose(fid);

%{
fid=fopen('neuinfo.dat','r');
ln=6;
nn=fscanf(fid,'%f',prod(lxs)*ln);
nn=reshape(nn,[lxs(1),lxs(2),lxs(3),ln]);
fclose(fid);
%}


%% REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if lxs(2) == 1 || lxs(3) == 1    %a 2D simulations was done
  %Jpar=squeeze(J1);
 % Jperp2=squeeze(J3);
  dat.ne = squeeze(dat.ns(:,:,:,lsp));
  %p=dat.ns(:,:,:,1)./dat.ns(:,:,:,lsp);
  %p=squeeze(p);
  %vi=squeeze(v1);
  %vi2=squeeze(v2);
%  vi3=permute(v3,[3,2,1]);
  dat.Ti = squeeze(sum(dat.ns(:,:,:,1:6) .* dat.Ts(:,:,:,1:6),4) ./ dat.ns(:,:,:,lsp));
  dat.Te = squeeze(dat.Ts(:,:,:,lsp));

%  [X3,X1]=meshgrid(x3,x1);
else    %full 3D run
%  Jpar=permute(J1(:,:,:),[3,2,1]);
%  Jperp2=permute(J2(:,:,:),[3,2,1]);
%  Jperp3=permute(J3(:,:,:),[3,2,1]);
%  ne=permute(dat.ns(:,:,:,lsp),[3,2,1]);
%  p=dat.ns(:,:,:,1)./dat.ns(:,:,:,lsp);
%  p=permute(p,[3,2,1]);
%  vi=permute(v1,[3,2,1]);
%  vi2=permute(v2,[3,2,1]);
%  vi3=permute(v3,[3,2,1]);
%  Ti=sum(dat.ns(:,:,:,1:6).*dat.Ts(:,:,:,1:6),4)./dat.ns(:,:,:,lsp);
%  Ti=permute(Ti,[3,2,1]);
%  Te=permute(dat.Ts(:,:,:,lsp),[3,2,1]);
%
%  [X2,X3,X1]=meshgrid(x2,x3,x1);

  %Jpar=J1(:,:,:);
  %Jperp2=J2(:,:,:);
  %Jperp3=J3(:,:,:);
  dat.ne = dat.ns(:,:,:,lsp);
  %p=dat.ns(:,:,:,1)./dat.ns(:,:,:,lsp);
  %p=p;
  %vi=v1;
  %vi2=v2;
  %vi3=v3;
  dat.Ti = sum(dat.ns(:,:,:,1:6) .* dat.Ts(:,:,:,1:6),4) ./ dat.ns(:,:,:,lsp);
  dat.Te = dat.Ts(:,:,:,lsp);
end

end % function
