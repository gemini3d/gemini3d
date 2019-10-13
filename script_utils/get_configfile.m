function filename = get_configfile(path)
%% get configuration file, checking first for config.nml and then config.ini

narginchk(1,1)

if is_file(path)
  filename = path;
elseif is_folder(path)
  names = {'config.nml', 'config.ini'};
  for s = names
    filename = [path, filesep, s{:}];
    if is_file(filename)
      break
    end
  end
else
  error(['could not find config file in ', path])
end

assert(is_file(filename), ['could not find config file in ', filename])

end % function