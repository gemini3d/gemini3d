function filename = get_configfile(path)
%% get configuration file

narginchk(1,1)

% necessary for Matlab
path = absolute_path(path);

if is_file(path)
  filename = path;
elseif is_folder(path)
  names = {'config.nml', 'inputs/config.nml', 'config.ini', 'inputs/config.ini'};
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
