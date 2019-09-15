function filename = get_configfile(path)

validateattr(path, {'char'}, {'vector'}, mfilename,'directory or file',1)

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