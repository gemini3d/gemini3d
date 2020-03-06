function dat = loadframe3Dcurvne(filename)

narginchk(1,1)

[~,~,ext] = fileparts(filename);

switch ext
  case '.dat'
  case {'.h5', '.nc'}, error('only raw is handled for now. Raise a Github Issue')
  otherwise, error(['unknown file type ',filename])
end

%% SIMULATION SIZE
lxs = simsize(direc);
%% SIMULATION RESULTS
assert(is_file(filename), [filename,' is not a file.'])

fid=fopen(filename,'r');
simdt(fid);

ns=fread(fid,prod(lxs),'real*8');
ns=reshape(ns, lxs);

fclose(fid);

%REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if (lxs(2) == 1)    %a 2D simulations was done
  dat.ne = squeeze(ns(:,:,:));
%  [X3,X1]=meshgrid(x3,x1);
else    %full 3D run
%  ne=permute(ns(:,:,:),[3,2,1]);
  dat.ne = ns;
%  [X2,X3,X1]=meshgrid(x2,x3,x1);
end

end % function
