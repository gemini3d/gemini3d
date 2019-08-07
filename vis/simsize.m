function lxs = simsize(direc)

validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)

fn = [direc,filesep,'inputs', filesep, 'simsize.dat'];
assert(exist(fn,'file')==2, [fn,' does not exist'])

fid = fopen(fn ,'r');
lxs = fread(fid,3,'integer*4');
lxs = lxs(:).';  % needed for concatenation
fclose(fid);

end