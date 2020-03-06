function dat = loadframe3Dcurv_hdf5(fn)

narginchk(1,1)
%% SIMULATION SIZE
lsp=7;
lxs = simsize(fn);
%% SIMULATION RESULTS
assert(is_file(fn), [fn,' is not a file.'])
dat.filename = fn;

if isoctave
  D = load(fn);

  dat.ns = D.nsall;
  dat.vs1 = D.vs1all;
  dat.Ts = D.Tsall;
  dat.J1 = D.J1all;
  dat.J2 = D.J2all;
  dat.J3 = D.J3all;
  dat.v2 = D.v2avgall;
  dat.v3 = D.v3avgall;
  dat.Phitop = D.Phiall;
else
  dat.ns = h5read(fn, '/nsall');
  dat.vs1 = h5read(fn, '/vs1all');
  dat.Ts = h5read(fn, '/Tsall');
  dat.J1 = h5read(fn, '/J1all');
  dat.J2 = h5read(fn, '/J2all');
  dat.J3 = h5read(fn, '/J3all');
  dat.v2 = h5read(fn, '/v2avgall');
  dat.v3 = h5read(fn, '/v3avgall');
  dat.Phitop = h5read(fn, '/Phiall');
end

dat.J1 = squeeze(dat.J1);
dat.J2 = squeeze(dat.J2);
dat.J3 = squeeze(dat.J3);
dat.v2 = squeeze(dat.v2);
dat.v3 = squeeze(dat.v3);

dat.v1 = squeeze(sum(dat.ns(:,:,:,1:6) .* dat.vs1(:,:,:,1:6), 4) ./ dat.ns(:,:,:,lsp));

%% REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if lxs(2) == 1 || lxs(3) == 1  %a 2D simulations was done
  dat.ne=squeeze(dat.ns(:,:,:,lsp));

  dat.Ti = squeeze(sum(dat.ns(:,:,:,1:6) .* dat.Ts(:,:,:,1:6),4) ./ dat.ns(:,:,:,lsp));
  dat.Te = squeeze(dat.Ts(:,:,:,lsp));

%  [X3,X1]=meshgrid(x3,x1);
else    %full 3D run
  dat.ne = dat.ns(:,:,:,lsp);

  dat.Ti = sum(dat.ns(:,:,:,1:6) .* dat.Ts(:,:,:,1:6),4) ./ dat.ns(:,:,:,lsp);
  dat.Te = dat.Ts(:,:,:,lsp);
end

end % function
