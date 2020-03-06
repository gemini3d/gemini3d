function makedir(path)
%% malformed paths can be "created" but are not accessible.
% This function workaround that bug in Matlab mkdir().
%
narginchk(1,1)

path = absolute_path(path);

if ~is_folder(path)
  mkdir(path);
end

assert(is_folder(path), [path, ' does not exist'])
end
