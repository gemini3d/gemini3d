function h5save(filename, varname, A, sizeA, dtype)

narginchk(3, 5)

if nargin < 4 || isempty(sizeA)
  if isvector(A)
    sizeA = length(A);
  else
    sizeA = size(A);
  end
end

% coerce if needed
if nargin >= 5 && ~isempty(dtype)
  switch dtype
    case {'float64', 'double'}
      if ~isa(A, 'double')
        A = double(A);
      end
    case {'float32', 'single'}
      if ~isa(A, 'single')
        A = single(A);
      end
    otherwise, error(['unknown dtype ',dtype])
  end
end

varnames = {};
if is_file(filename)
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
  if ~ismatrix(A)
    % enable Gzip compression--remember Matlab's dim order is flipped from
    % C / Python
    switch ndims(A)
      case 4, chunksize = [sizeA(1), sizeA(2), 1, sizeA(4)];
      case 3, chunksize = [sizeA(1), sizeA(2), 1];
      otherwise, error('FIXME: bigger than 4 dims')
    end
    h5create(filename, varname, sizeA, 'DataType', class(A), ...
      'Deflate', 1, 'Fletcher32', true, 'Shuffle', true, 'ChunkSize', chunksize)
  else
    h5create(filename, varname, sizeA, 'DataType', class(A))
  end
  h5write(filename, varname, A)
end

end
