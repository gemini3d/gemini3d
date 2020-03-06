function p = read_config(path)
% reads simulation configuration into struct
narginchk(1,1)

filename = get_configfile(path);

[~,~,ext] = fileparts(filename);

switch ext
  case '.nml', p = read_nml(filename);
  case '.ini', p = read_ini(filename);
  otherwise, error(['not sure how to read config file ', filename])
end

%% deduce data file format from simsize format
[~,~,ext] = fileparts(p.indat_size);
p.format = ext(2:end);

end % function