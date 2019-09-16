%!assert(isstruct(read_nml_group))
function params = read_nml_group(filename, group)

narginchk(2,2)
assert(is_file(filename), ['config file ', filename, ' not found.'])
validateattr(group, {'char'}, {'vector'}, mfilename, 'nml group name', 2)

params = struct();

fid=fopen(filename);
while ~feof(fid)
  line = fgetl(fid);
  mgrp = regexp(line, ['^\s*&', group], 'match');
  if ~isempty(mgrp)
     break
  end
end

assert(~isempty(mgrp), ['did not find group ', group, ' in ', filename])

while ~feof(fid)
  line = fgetl(fid);
  mend = regexp(line, '^\s*/\s*$', 'match');
  if ~isempty(mend)
     break
  end

  [k, v] = strtok(line, '=');
  v = strtok(v(2:end), '!');  % discard comments
  % is it floats of up to length 3?
  vals = sscanf(v, '%f,%f,%f');
  if isempty(vals)  % must be a string
    vals = v;
  end

  params.(strtrim(k)) = vals;
end
fclose(fid);

assert(~isempty(mend), ['did not read end of group ', group, ' in ', filename])

end % function