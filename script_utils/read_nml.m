function params = read_nml(filename)

narginchk(1,1)

base = read_nml_group(filename, 'base');

params = base;

end % function