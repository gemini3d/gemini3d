function writedata(dmy,time,ns,vsx1,Ts,outdir,outID)

%--------------------------------------------------------
%-----WRITE STATE VARIABLE DATA TO BE USED AS INITIAL
%-----CONDITIONS FOR ANOTHER SIMULATION.  NOTE THAT WE
%-----DO NOT HERE OUTPUT ANY OF THE ELECTRODYNAMIC
%-----VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
%-----UP IN THE FORTRAN CODE.
%-----
%-----INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
%-----(I.E. THEY SHOULD NOT INCLUDE GHOST CELLS
%--------------------------------------------------------


  %% CHECK THAT THE DIRECTORY EXISTS - IF NOT CREATE IT AND GIVE THE USER AN INDICATION OF WHAT WAS DONE
  if (~(exist(outdir,'dir')==7))
    mkdir(outdir);
    disp(['Note that directory:  ',outdir,' needed to be created for sim. data file']);
  end


  %% OPEN THE FILE
  filename=[outdir,'/',outID,'_ICs.dat'];
  fid=fopen(filename,'w');

  
  %% WRITE THE DATA
  fwrite(fid,dmy,'real*8');
  fwrite(fid,time,'real*8');
  fwrite(fid,ns,'real*8');
  fwrite(fid,vsx1,'real*8');
  fwrite(fid,Ts,'real*8');

  fclose(fid);
  
end
