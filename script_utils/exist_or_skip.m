function exist_or_skip(filename, dirfile)
% if paths not exist, exit code 77 as GNU standard skip test indicator
narginchk(2,2)
validateattr(dirfile, {'char'}, {'vector'})

if strcmp(dirfile, 'file')
    if ~is_file(filename)
      if isinteractive
        error([filename, ' is not a file.'])
      else
        fprintf(2, [filename, ' is not a file.\n'])
        exit(77)
      endif
    end
elseif strcmp(dirfile, 'dir')
    if ~is_folder(filename)
      if isinteractive
        error([filename, ' is not a directory.'])
      else
        fprintf(2, [filename, ' is not a directory.\n'])
        exit(77)
      endif
    end

end

end % function
