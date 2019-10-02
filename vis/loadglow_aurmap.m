function cAur = loadglow_aurmap(filename, lx2, lx3, lwave)
narginchk(4,4)
assert(is_file(filename), [filename, " is not a filename"])
validateattributes(lx2, {'numeric'}, {'scalar', 'integer', 'positive'})
validateattributes(lx3, {'numeric'}, {'scalar', 'integer', 'positive'})
validateattributes(lwave, {'numeric'}, {'scalar', 'integer', 'positive'})

fid=fopen(filename,'r');
cAur = fread(fid,lx2*lx3*lwave,'real*8');

cAur = reshape(cAur,[lx2,lx3,lwave]);
%cAur = reshape(cAur,[lwave,lx2,lx3]);
%cAur=permute(cAur,[2,3,1]);

fclose(fid);

end