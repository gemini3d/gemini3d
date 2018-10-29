function dat = read4D(fid, lsp, lxs)

dat = fread(fid, prod(lxs)*lsp,'real*8');

dat = reshape(dat, [lxs,lsp]);

end
