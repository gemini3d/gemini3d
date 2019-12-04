function abspath = absolute_path(path)
% path need not exist, but absolute path is returned

% have to expand ~ first
path = expanduser(path);

if isoctave
  abspath = make_absolute_filename(path);
else
  abspath = char(java.io.File(path).getCanonicalPath());
end

end