%!assert(ischar(absolute_path('~')))
function abspath = absolute_path(path)
% path need not exist, but absolute path is returned
%
% NOTE: Octave is weaker at this, especially if /foo/bar/../baz and "bar"
% doesn't exist, it may just return the input unmodified.

narginchk(1,1)
% have to expand ~ first
path = expanduser(path);

if isoctave
  abspath = make_absolute_filename(path);
else
  if is_relative_path(path)
    % TODO: handle .file or .path/.file
    path = [pwd, filesep, path];
  end
  abspath = char(java.io.File(path).getCanonicalPath());
end

% debugging
% disp([path, ' => ',abspath])

end % function
