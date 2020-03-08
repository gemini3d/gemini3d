function [ymd,UTsec,tdur,dtout,flagoutput,mloc,activ,indat_size,indat_grid,indat_file] = readconfig(path)

narginchk(1,1)

params = read_config(path);


%% FIXME: make rest of programs use struct instead of individual variables
ymd = params.ymd;
UTsec = params.UTsec0;
tdur = params.tdur;
dtout = params.dtout;
flagoutput = params.flagoutput;
activ = params.activ;
try
  mloc = params.mloc;
catch
  mloc = [];
end
indat_size=params.indat_size;
indat_grid=params.indat_grid;
indat_file=params.indat_file;
end % function