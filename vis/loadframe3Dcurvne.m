function ne = loadframe3Dcurvne(direc, filename)

%SIMULATION SIZE
fid=fopen([direc,'/inputs/simsize.dat'],'r');
lxs=fread(fid,3,'integer*4');
lxs=lxs(:)';
fclose(fid);


%SIMULATION GRID FILE (NOTE THAT THIS IS NOT THE ENTIRE THING - THAT NEEDS
%TO BE DONE WITH READGRID.M.  WE NEED THIS HERE TO DO MESHGRIDS
lsp=7;

%THESE ARE THE SIMULATIONS RESULTS
fid=fopen([direc,filename],'r');
%t=fread(fid,1,'real*8');     %this is for just UT seconds
simdate=zeros(1,6);    %datevec-style array
simdatetmp=fread(fid,4,'real*8');
simdate(1:4)=simdatetmp;

ns=fread(fid,prod(lxs),'real*8');
ns=reshape(ns,[lxs]);
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
