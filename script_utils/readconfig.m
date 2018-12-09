function [ymd,UTsec,tdur,dtout,flagoutput,mloc] = readconfig(filename)

  validateattr(filename, {'char'}, {'vector'}, mfilename, 'configuration filename', 1)
  
  assert(exist(filename,'file')==2, [filename,' does not exist'])
    
  fid=fopen(filename);

  %DATE
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  [day,remainder]=strtok(datatrim,',');
  [month,remainder]=strtok(remainder,',');
  [year,remainder]=strtok(remainder,',');    %should not find delimiter..
  ymd=[str2num(year),str2num(month),str2num(day)];

  %UT seconds
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  UTsec=str2num(datatrim);

  %Sim duration
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  tdur=str2num(datatrim);  
 
  %Output dt
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  dtout=str2num(datatrim);

  %Strip out the junk
%  data=fgetl(fid);
  data=fgetl(fid);
%  data=fgetl(fid);
%  data=fgetl(fid);

  %CFL number
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  tcfl=str2num(datatrim); 

  %Exospheric temp.
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  Teinf=str2num(datatrim);

  %Strip out the junk
%  data=fgetl(fid);

  %Type of potential solution
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  flagpot=str2num(datatrim);

  %Grid periodic flag
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  flagperiodic=str2num(datatrim);

  %OUTPUT type flag
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  flagoutput=str2num(datatrim);

  %capacitance flag
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  flagcap=str2num(datatrim);

  %Strip out the junk
  data=fgetl(fid);
  data=fgetl(fid);
  data=fgetl(fid);

  %Neutral perturbation info
  data=fgetl(fid);
  [datatrim,remainder]=strtok(data,' ');
  flagdneu=str2num(datatrim);

  %Type of neutral interpolation done (not used)
  data=fgetl(fid);

  if (flagdneu)
    data=fgetl(fid);
    [datatrim,remainder]=strtok(data,' ');
    [mlat,remainder]=strtok(datatrim,',');
    [mlon,remainder]=strtok(remainder,',');
    mloc=[str2num(mlat),str2num(mlon)];
  else
    mloc=[];
  end
  
  fclose(fid);

  %There's more in the input file, but we don't use it in the data processing
end
