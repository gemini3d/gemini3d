function [ymd,UTsec,tdur,dtout,flagoutput,mloc,activ] = readconfig(path)

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
end % function