function [x1, x2, x3] = simaxes(direc)

narginchk(1,1)
fn = [direc,'/inputs/simgrid.dat'];
assert(is_file(fn), [fn, ' is not a file.'])

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