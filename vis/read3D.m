function dat = read3D(fid, lxs)
narginchk(2,2)
dat = fread(fid, lxs(1)*lxs(2)*lxs(3), 'real*8');
dat = reshape(dat, [lxs(1),lxs(2),lxs(3)]);
end
