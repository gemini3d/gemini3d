function dat = read4D(fid, lsp, lxs)
narginchk(3,3)
dat = fread(fid, prod(lxs)*lsp,'real*8');

dat = reshape(dat, [lxs,lsp]);

end
