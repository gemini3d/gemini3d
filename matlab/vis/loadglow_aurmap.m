function cAur = loadglow_aurmap(filename, lx2, lx3, lwave)
%% loads simulated auroral emissions
narginchk(4,4)
assert(is_file(filename), [filename, " is not a filename"])
validateattributes(lx2, {'numeric'}, {'scalar', 'integer', 'positive'})
validateattributes(lx3, {'numeric'}, {'scalar', 'integer', 'positive'})
validateattributes(lwave, {'numeric'}, {'scalar', 'integer', 'positive'})

[~,~,ext] = fileparts(filename);

switch ext
  case '.dat', cAur = loadglow_aurmap_raw(filename, lx2, lx3, lwave);
  case '.h5', cAur = loadglow_aurmap_hdf5(filename);
  otherwise, error(['unknown file type ',filename])
end

end


function cAur = loadglow_aurmap_raw(filename, lx2, lx3, lwave)
narginchk(4,4)

fid=fopen(filename,'r');
cAur = fread(fid,lx2*lx3*lwave,'real*8');
fclose(fid);

cAur = reshape(cAur,[lx2,lx3,lwave]);
%cAur = reshape(cAur,[lwave,lx2,lx3]);
%cAur=permute(cAur,[2,3,1]);
end % function


function cAur = loadglow_aurmap_hdf5(filename)
narginchk(1,1)

if isoctave
  D = load(filename);
  cAur = D.aurora.iverall;
else
  cAur = h5read(filename, '/aurora/iverout');
end

cAur = squeeze(cAur);

end % function