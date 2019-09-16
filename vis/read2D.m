function dat = read2D(fid, lxs)
narginchk(2,2)
dat = fread(fid, lxs(2)*lxs(3), 'real*8');
dat = reshape(dat, [lxs(2),lxs(3)]);
end
