function [ne,v1,Ti,Te,J1,v2,v3,J2,J3,ns,vs1,Ts,Phitop] = loadframe3Dcurv_hdf5(fn)

narginchk(1,1)
%% SIMULATION SIZE
lsp=7;
lxs = simsize(fn);
%% SIMULATION RESULTS
assert(is_file(fn), [fn,' is not a file.'])

if isoctave
  D = load(fn);

  ns = D.nsall;
  vs1 = D.vs1all;
  Ts = D.Tsall;
  J1 = D.J1all;
  J2 = D.J2all;
  J3 = D.J3all;
  v2 = D.v2avgall;
  v3 = D.v3avgall;
  Phitop = D.Phiall;
else
  ns = h5read(fn, '/nsall');
  vs1 = h5read(fn, '/vs1all');
  Ts = h5read(fn, '/Tsall');
  J1 = h5read(fn, '/J1all');
  J2 = h5read(fn, '/J2all');
  J3 = h5read(fn, '/J3all');
  v2 = h5read(fn, '/v2avgall');
  v3 = h5read(fn, '/v3avgall');
  Phitop = h5read(fn, '/Phiall');
end

v1 = sum(ns(:,:,:,1:6) .* vs1(:,:,:,1:6), 4) ./ ns(:,:,:,lsp);

%% REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if lxs(2) == 1   %a 2D simulations was done
  ne=squeeze(ns(:,:,:,lsp));

  Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
  Ti=squeeze(Ti);
  Te=squeeze(Ts(:,:,:,lsp));

%  [X3,X1]=meshgrid(x3,x1);
else    %full 3D run
  ne=ns(:,:,:,lsp);

  Ti=sum(ns(:,:,:,1:6).*Ts(:,:,:,1:6),4)./ns(:,:,:,lsp);
  Te=Ts(:,:,:,lsp);
end

end % function
