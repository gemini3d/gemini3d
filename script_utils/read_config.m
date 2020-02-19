function params = read_config(path)
% reads simulation configuration into struct
narginchk(1,1)

filename = get_configfile(path);

[~,~,ext] = fileparts(filename);

switch ext
  case '.ini', params = read_ini(filename);
  case '.nml', params = read_nml(filename);
  otherwise, error(['not sure how to read config file ', filename])
end

end % function