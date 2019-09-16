function lxs = simsize(direc)

narginchk(1,1)
assert(is_folder(direc), [direc, ' is not a directory.'])

fn = [direc,filesep, 'inputs', filesep, 'simsize.dat'];
assert(is_file(fn), [fn,' is not a file.'])

fid = fopen(fn, 'r');
lxs = fread(fid, 3, 'integer*4');
lxs = lxs(:).';  % needed for concatenation
fclose(fid);

end