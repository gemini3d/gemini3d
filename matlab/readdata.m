function [t,ns,Ts,vs1,J1,J2,J3,v2,v3,Phitop]=readdata(lxs,filename)

narginchk(2,2)
validateattr(lxs, {'numeric'}, {'vector', 'numel', 3, 'positive'}, mfilename, 'grid dimensions', 1)

filename = absolute_path(filename);
assert(is_file(filename), [filename, 'is not found.'])
%% READ DATA FROM AN OUTPUT FILE WRITTEN BY FORTRAN CODE

fid=fopen(filename,'r');
t=fread(fid,1,'real*8');

ns=fread(fid,prod(lxs)*lsp,'real*8');
ns=reshape(ns,[lxs,lsp]);
disp('READDATA --> Loaded densities...')

vs1=fread(fid,prod(lxs)*lsp,'real*8');
vs1=reshape(vs1,[lxs,lsp]);
v1=sum(ns(:,:,:,1:6).*vs1(:,:,:,1:6),4)./ns(:,:,:,lsp);
disp('READDATA --> Loaded parallel velocities...')

Ts=fread(fid,prod(lxs)*lsp,'real*8');
Ts=reshape(Ts,[lxs,lsp]);
disp('READDATA --> Loaded temperatures...'))

if (~feof(fid))   %some files may not have electrodynamic info
    J1=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
    J1=reshape(J1,[lxs(1),lxs(2),lxs(3)]);

    J2=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
    J2=reshape(J2,[lxs(1),lxs(2),lxs(3)]);

    J3=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
    J3=reshape(J3,[lxs(1),lxs(2),lxs(3)]);
    disp('READDATA --> Loaded current density...')

    v2=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
    v2=reshape(v2,[lxs(1),lxs(2),lxs(3)]);

    v3=fread(fid,lxs(1)*lxs(2)*lxs(3),'real*8');
    v3=reshape(v3,[lxs(1),lxs(2),lxs(3)]);
    disp('READDATA --> Loaded perpendicular drifts...')

    Phitop=fread(fid,lxs(2)*lxs(3),'real*8');
    Phitop=reshape(Phitop,[lxs(2),lxs(3)]);
    disp('READDATA --> Loaded topside potential pattern...')
else
    J1=[];
    J2=[];
    J3=[];
    v2=[];
    v3=[];
    Phitop=[];
    disp('READDATA --> Skipping electrodynamic parameters...')
end

fclose(fid);


%REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if (lxs(2) == 1)    %a 2D simulations was done
    Jpar=squeeze(J1);
    Jperp2=squeeze(J3);
    ne=squeeze(ns(:,:,:,lsp));
    p=ns(:,:,:,1)./ns(:,:,:,lsp);
    p=squeeze(p);
    vi=squeeze(v1);
    vi2=squeeze(v2);
    %  vi3=permute(v3,[3,2,1]);
    Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
    Ti=squeeze(Ti);
    Te=squeeze(Ts(:,:,:,lsp));

    [X3,X1]=meshgrid(x3,x1);
else    %full 3D run
    Jpar=permute(J1(:,:,:),[3,2,1]);
    Jperp2=permute(J2(:,:,:),[3,2,1]);
    Jperp3=permute(J3(:,:,:),[3,2,1]);
    ne=permute(ns(:,:,:,lsp),[3,2,1]);
    p=ns(:,:,:,1)./ns(:,:,:,lsp);
    p=permute(p,[3,2,1]);
    vi=permute(v1,[3,2,1]);
    vi2=permute(v2,[3,2,1]);
    vi3=permute(v3,[3,2,1]);
    Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
    Ti=permute(Ti,[3,2,1]);
    Te=permute(Ts(:,:,:,lsp),[3,2,1]);

    [X2,X3,X1]=meshgrid(x2,x3,x1);
end

end
