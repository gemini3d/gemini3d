function [ymd,UTsec,tdur,dtout,flagoutput,mloc,activ] = readconfig(path)

narginchk(1,1)
validateattr(path, {'char'}, {'vector'}, mfilename, 'configuration path/filename', 1)

assert(is_folder(path) || is_file(path), ['configuration path ', path,' does not exist'])

filename = get_configfile(path);

[~,~,ext] = fileparts(filename);

switch ext
    case '.ini', params = read_ini(filename);
    case '.nml', params = read_nml(filename);
    otherwise, error(['not sure how to read config file ', filename])
end

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