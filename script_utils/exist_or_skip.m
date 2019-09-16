function exist_or_skip(filename)
% if paths not exist, exit code 77 as GNU standard skip test indicator

if ~is_file(filename)
  if isinteractive
    error([filename, ' not found'])
  else
    fprintf(2, [filename, ' not found\n'])
    exit(77)
  endif
end

end % function
