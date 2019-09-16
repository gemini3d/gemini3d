function ne = loadframe3Dcurvne(direc, filename)

narginchk(2,2)
%% SIMULATION SIZE
lxs = simsize(direc);
%% SIMULATION RESULTS
fsimres = [direc,filesep,filename];
assert(is_file(fsimres), [fsimres,' is not a file.'])

fid=fopen(fsimres,'r');
simdt(fid);

ns=fread(fid,prod(lxs),'real*8');
ns=reshape(ns, lxs);
fprintf('Loaded densities...\n');

fclose(fid);

%REORGANIZE ACCORDING TO MATLABS CONCEPT OF A 2D or 3D DATA SET
if (lxs(2) == 1)    %a 2D simulations was done
  ne=squeeze(ns(:,:,:));
%  [X3,X1]=meshgrid(x3,x1);
else    %full 3D run
%  ne=permute(ns(:,:,:),[3,2,1]);
  ne=ns;
%  [X2,X3,X1]=meshgrid(x2,x3,x1);
end

end % function
