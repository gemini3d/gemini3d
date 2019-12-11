%!assert(ischar(absolute_path('~')))
function abspath = absolute_path(path)
% path need not exist, but absolute path is returned

% have to expand ~ first
path = expanduser(path);

if isoctave
  abspath = make_absolute_filename(path);
else
  if path(1) == '.'
    % TODO: handle .file or .path/.file
    % this is necessary when run from arbitrary directories.
    path = [pwd, filesep, path];
  end
  abspath = char(java.io.File(path).getCanonicalPath());
end

end % function
