%!assert(ischar(absolute_path('~')))
function abspath = absolute_path(path)
% path need not exist, but absolute path is returned

narginchk(1,1)
% have to expand ~ first
path = expanduser(path);

if isoctave
  abspath = make_absolute_filename(path);
else
  if path(1) == '.' || ...
     ~contains(path, filesep) || ...
     (isletter(path(1)) && ~strcmp(path(2), ':'))
    % TODO: handle .file or .path/.file
    path = [pwd, filesep, path];
  end
  abspath = char(java.io.File(path).getCanonicalPath());
end

end % function
