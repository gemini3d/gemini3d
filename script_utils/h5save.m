function h5save(filename, varname, A, sizeA)

narginchk(3, 4)

if nargin < 4
  if isvector(A)
    sizeA = length(A);
  else
    sizeA = size(A);
  end
end

varnames = {};
if isfile(filename)
  varnames = extractfield(h5info(filename).Datasets, 'Name');
end

if any(strcmp(varname, varnames) | strcmp(varname(2:end), varnames))
  % existing variable
  diskshape = h5info(filename, varname).Dataspace.Size;
  if length(diskshape) >= 2
    if diskshape(1) == 1 % isrow
      start = ones(ndims(A),1);
    elseif diskshape(2) == 1 % iscolumn
      start = ones(1,ndims(A));
    else
      start = ones(1,ndims(A));
    end
  else
    start = 1;
  end

  if all(diskshape == sizeA)
    h5write(filename, varname, A, start, sizeA)
  elseif all(diskshape == fliplr(sizeA))
    h5write(filename, varname, A.', start, fliplr(sizeA))
  else
    error(['shape of ',varname,': ',int2str(sizeA),' does not match existing HDF5 shape: ', int2str(diskshape)])
  end
else % new variable
  h5create(filename, varname, sizeA, 'DataType', class(A))
  h5write(filename, varname, A)
end

end
