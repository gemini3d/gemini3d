function ne = loadframe3Dcurvne(direc, filename)

validateattr(direc, {'char'}, {'vector'}, mfilename, 'data directory', 1)
validateattr(direc, {'char'}, {'vector'}, mfilename, 'data filename', 2)
%% SIMULATION SIZE
lxs = simsize(direc);
disp(['sim grid dimensions: ',num2str(lxs)])
%% SIMULATION RESULTS
fsimres = [direc,filesep,filename];
assert(exist(fsimres,'file')==2, [fsimres,' does not exist'])

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
