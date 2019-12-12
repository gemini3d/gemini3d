function h5save(filename, varname, A, sizeA)
narginchk(3, 4)
if nargin < 4, sizeA = size(A); end

varnames = {};
if isfile(filename)
  varnames = extractfield(h5info(filename).Datasets, 'Name');
end

if any(strcmp(varname, varnames))
  % FIXME: existing variable
  h5write(filename, varname, A, 1, length(A))
else % new variable
  h5create(filename, varname, sizeA)
  h5write(filename, varname, A)
end

end
