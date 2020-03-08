function params = read_ini(filename)
% read config*.ini files. Only for legacy use, new simulations should used
% config.nml Namelist format.
narginchk(1,1)

assert(is_file(filename), ['config file ', filename, ' not found.'])

fid=fopen(filename);

%DATE
datatrim = strtok(fgetl(fid),' ');
[day,remainder]=strtok(datatrim,',');
[month,remainder]=strtok(remainder,',');
year =strtok(remainder,',');    %should not find delimiter..
params.ymd=[str2double(year),str2double(month),str2double(day)];

%UT seconds
datatrim =strtok(fgetl(fid),' ');
params.UTsec0=str2double(datatrim);

%Sim duration
datatrim =strtok(fgetl(fid),' ');
params.tdur=str2double(datatrim);

%Output dt
datatrim =strtok(fgetl(fid),' ');
params.dtout=str2double(datatrim);

%f10.7 and geomag indices (used in msis)
datatrim =strtok(fgetl(fid),' ');
[f107a,remainder]=strtok(datatrim,',');
[f107,remainder]=strtok(remainder,',');
ap =strtok(remainder,',');    %should not find delimiter..
params.activ=[str2double(f107a),str2double(f107),str2double(ap)];

%CFL number
datatrim =strtok(fgetl(fid),' ');
params.tcfl=str2double(datatrim);

%Exospheric temp.
datatrim =strtok(fgetl(fid),' ');
params.Teinf=str2double(datatrim);

%Strip out the junk
%  data=fgetl(fid);

%Type of potential solution
datatrim = strtok(fgetl(fid),' ');
params.flagpot=str2double(datatrim);

%Grid periodic flag
datatrim =strtok(fgetl(fid),' ');
params.flagperiodic=str2double(datatrim);

%OUTPUT type flag
datatrim = strtok(fgetl(fid),' ');
params.flagoutput=str2double(datatrim);

%capacitance flag
datatrim = strtok(fgetl(fid),' ');
params.flagcap=str2double(datatrim);

%names of the input files (to be) used in the simulation
datatrim=fgetl(fid);
params.indat_size=datatrim;
datatrim=fgetl(fid);
params.indat_grid=datatrim;
datatrim=fgetl(fid);
params.indat_file=datatrim;

%Neutral perturbation info
datatrim = strtok(fgetl(fid),' ');
params.flagdneu = str2double(datatrim);

%Type of neutral interpolation done (not used)
fgetl(fid);

params.mloc=[];
if params.flagdneu
  datatrim = strtok(fgetl(fid),' ');
  [mlat,remainder]=strtok(datatrim,',');
  mlon = strtok(remainder,',');
  % NOTE: str2num necessary here in case input formatted like 1.3d0 which
  % str2double returns NaN for while str2num works
  params.mloc=[str2num(mlat), str2num(mlon)]; %#ok<ST2NM>
end

fclose(fid);

%There's more in the input file, but we don't use it in the data processing
end % function
