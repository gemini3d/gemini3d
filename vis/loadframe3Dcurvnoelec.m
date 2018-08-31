function [ne,v1,Ti,Te,ns,Ts,vs1,simdate] = loadframe3Dcurvnoelec(direc, filename)

%READ IN SIMULATION DATA WITH NO ELECTRODYNAMIC PARAMS SUCH AS FROM AN
%INPUT FILE

%SIMULATION SIZE
fid=fopen([direc,'simsize.dat'],'r');
lxs=fread(fid,3,'integer*4');
lxs=lxs(:)';
fclose(fid);


%SIMULATION GRID FILE (NOTE THAT THIS IS NOT THE ENTIRE THING - THAT NEEDS
%TO BE DONE WITH READGRID.M.  WE NEED THIS HERE TO DO MESHGRIDS
lsp=7;
fid=fopen([direc,'simgrid.dat'],'r');
x1=fread(fid,lxs(1),'real*8');
x1=x1(:)';
x2=fread(fid,lxs(2),'real*8');
x2=x2(:)';
x3=fread(fid,lxs(3),'real*8');
x3=x3(:)';
fclose(fid);


%THESE ARE THE SIMULATIONS RESULTS
fid=fopen([direc,filename],'r');
%t=fread(fid,1,'real*8');     %this is for just UT seconds
simdate=zeros(1,6);    %datevec-style array
simdatetmp=fread(fid,4,'real*8');
simdate(1:4)=simdatetmp;

ns=fread(fid,prod(lxs)*lsp,'real*8');
ns=reshape(ns,[lxs,lsp]);
fprintf('Loaded densities...\n');

vs1=fread(fid,prod(lxs)*lsp,'real*8');
vs1=reshape(vs1,[lxs,lsp]);
v1=sum(ns(:,:,:,1:6).*vs1(:,:,:,1:6),4)./ns(:,:,:,lsp);
fprintf('Loaded parallel velocities...\n');

Ts=fread(fid,prod(lxs)*lsp,'real*8');
Ts=reshape(Ts,[lxs,lsp]);
fprintf('Loaded temperatures...\n');

% J1=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
% J1=reshape(J1,[lxs(1),lxs(2),lxs(3)]);
% 
% J2=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
% J2=reshape(J2,[lxs(1),lxs(2),lxs(3)]);
% 
% J3=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
% J3=reshape(J3,[lxs(1),lxs(2),lxs(3)]);
% fprintf('Loaded current density...\n');
% 
% v2=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
% v2=reshape(v2,[lxs(1),lxs(2),lxs(3)]);
% 
% v3=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
% v3=reshape(v3,[lxs(1),lxs(2),lxs(3)]);
% fprintf('Loaded perpendicular drifts...\n');
% 
% Phitop=fread(fid,lxs(2)*lxs(3),'real*8');
% Phitop=reshape(Phitop,[lxs(2),lxs(3)]);
% fprintf('Loaded topside potential pattern...\n');

fclose(fid);

%{
fid=fopen('neuinfo.dat','r');
ln=6;
nn=fscanf(fid,'%f',prod(lxs)*ln);
nn=reshape(nn,[lxs(1),lxs(2),lxs(3),ln]);
fclose(fid);
%}


%REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if (lxs(2) == 1)    %a 2D simulations was done
%   Jpar=squeeze(J1);
%   Jperp2=squeeze(J3);
  ne=squeeze(ns(:,:,:,lsp));
  p=ns(:,:,:,1)./ns(:,:,:,lsp);
  p=squeeze(p);
%   vi=squeeze(v1);
%   vi2=squeeze(v2);
% %  vi3=permute(v3,[3,2,1]);
  Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
  Ti=squeeze(Ti);
  Te=squeeze(Ts(:,:,:,lsp));

%  [X3,X1]=meshgrid(x3,x1);
else    %full 3D run 
%  Jpar=permute(J1(:,:,:),[3,2,1]);
%  Jperp2=permute(J2(:,:,:),[3,2,1]);
%  Jperp3=permute(J3(:,:,:),[3,2,1]);
%  ne=permute(ns(:,:,:,lsp),[3,2,1]);
%  p=ns(:,:,:,1)./ns(:,:,:,lsp);
%  p=permute(p,[3,2,1]);
%  vi=permute(v1,[3,2,1]);
%  vi2=permute(v2,[3,2,1]);
%  vi3=permute(v3,[3,2,1]);
%  Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
%  Ti=permute(Ti,[3,2,1]);
%  Te=permute(Ts(:,:,:,lsp),[3,2,1]);
%
%  [X2,X3,X1]=meshgrid(x2,x3,x1);

%   Jpar=J1(:,:,:);
%   Jperp2=J2(:,:,:);
%   Jperp3=J3(:,:,:);
  ne=ns(:,:,:,lsp);
  p=ns(:,:,:,1)./ns(:,:,:,lsp);
  p=p;
  vi=v1;
%   vi2=v2;
%   vi3=v3;
  Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
  Te=Ts(:,:,:,lsp);
end
fprintf('Data organized...\n')



