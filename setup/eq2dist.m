function [nsi,vs1i,Tsi,xgin,ns,vs1,Ts]=eq2dist(eqdir,simID,xg)

  %READ IN THE SIMULATION INFORMATION
  [ymd0,UTsec0,tdur,dtout,flagoutput,mloc] = readconfig(eqdir, filesep, 'inputs');
  xgin = readgrid([eqdir,filesep, 'inputs/']);


  %FIND THE DATE OF THE END FRAME OF THE SIMULATION (PRESUMABLY THIS WILL BE THE STARTING POITN FOR ANOTEHR)
  [ymdend,UTsecend]=dateinc(tdur,ymd0,UTsec0);


  %LOAD THE FRAME
  [ne,mlatsrc,mlonsrc,xgin,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts]= ...
      loadframe(eqdir,ymdend,UTsecend,flagoutput,mloc,xgin);


  %DO THE INTERPOLATION
  [nsi,vs1i,Tsi]=model_resample(xgin,ns,vs1,Ts,xg);


  %THROW A WARNING IF THE INTERPOLATION WENT WEIRD...
  isfin_ns=min(isfinite(nsi(:)));
  isfin_vs1=min(isfinite(vs1i(:)));
  isfin_Ts=min(isfinite(Tsi(:)));
  if (~isfin_ns)
    warning('Detected a non-finite value in density...  Please fix before running your simulation...');
  end
  if (~isfin_vs1)
    warning('Detected a non-finite value in drift...  Please fix before running your simulation...');
  end
  if (~isfin_Ts)
    warning('Detected a non-finite value in temperature...  Please fix before running your simulation...');
  end


  %WRITE OUT THE GRID
  basedir=[eqdir,'../input/'];
  outdir=[basedir,simID];
  writegrid(xg,outdir);
  dmy=[ymdend(3),ymdend(2),ymdend(1)];
  writedata(dmy,UTsecend,nsi,vs1i,Tsi,outdir,simID);    %this uses the SIMID as both the output directory and as the filename tag...

end %function eq2dist
