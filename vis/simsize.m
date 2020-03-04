function lxs = simsize(path)

narginchk(1,1)

fn = [];
if is_file(path)
  [~, stem, ext] = fileparts(path);
  if strcmp(stem, 'simsize')
    fn = path;
  else
    fn = [fileparts(path), '/inputs/simsize', ext];
  end
elseif is_folder(path)
  for ext = {'.h5', '.dat'}
    fn = [path, '/inputs/simsize',ext{:}];
    if is_file(fn)
      break
    end
  end
end

assert(is_file(fn), [fn,' is not a file.'])
[~,~,ext] = fileparts(fn);

switch ext
  case '.h5'
    if isoctave
      d = load(fn);
      if isfield(d, 'lxs')
        lxs = d.lxs;
      elseif isfield(d, 'lx')
        lxs = d.lx;
      elseif isfield(d, 'lx1')
        % octave bug: octave_base_value::int32_scalar_value(): wrong type argument 'int32 matrix'
        lxs = [d.lx1; d.lx2; d.lx3];
      else
        error(['did not find lxs, lx, lx1 in ', fn])
      end
    else
      varnames = extractfield(h5info(fn).Datasets, 'Name');

      if any(strcmp('lxs', varnames))
        lxs = h5read(fn, '/lxs');
      elseif any(strcmp('lx', varnames))
        lxs = h5read(fn, '/lx');
      elseif any(strcmp('lx1', varnames))
        lxs = [h5read(fn, '/lx1'), h5read(fn, '/lx2'), h5read(fn, '/lx3')];
      else
        error(['did not find lxs, lx, lx1 in ', fn])
      end
    end
  case '.nc'
    varnames = extractfield(ncinfo(fn).Variables, 'Name');

    if any(strcmp('lxs', varnames))
      lxs = ncread(fn, '/lxs');
    elseif any(strcmp('lx', varnames))
        lxs = h5read(fn, '/lx');
    elseif any(strcmp('lx1', varnames))
      lxs = [ncread(fn, '/lx1'), ncread(fn, '/lx2'), ncread(fn, '/lx3')];
    end
  case '.dat'
    fid = fopen(fn, 'r');
    lxs = fread(fid, 3, 'integer*4');
    fclose(fid);
  otherwise, error(['unknown simsize file type ',fn])
end

lxs = lxs(:).';  % needed for concatenation

end % function
