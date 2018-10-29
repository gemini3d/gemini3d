function [x1, x2, x3] = simaxes(direc)

validateattributes(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)

fn = [direc,filesep,'inputs', filesep, 'simgrid.dat'];
assert(exist(fn,'file')==2, [fn,' does not exist'])

lxs = simsize(direc);

fid=fopen(fn, 'r');

x1=fread(fid,lxs(1),'real*8');
x1=x1(:)';

x2=fread(fid,lxs(2),'real*8');
x2=x2(:)';

x3=fread(fid,lxs(3),'real*8');
x3=x3(:)';

fclose(fid);

end