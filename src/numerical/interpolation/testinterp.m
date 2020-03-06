function testinterp(filename)

narginchk(1,1)

addpath([fileparts(mfilename('fullpath')), '/../../../matlab'])

exist_or_skip(filename, 'file')

if isoctave
  h = load(filename);
  lx1 = h.lx1;
  lx2 = h.lx2;
  x1 = h.x1;
  x2 = h.x2;
  f = h.f;
else
  lx1 = h5read(filename, '/lx1');
  lx2 = h5read(filename, '/lx2');
  x1 = h5read(filename, '/x1');
  x2 = h5read(filename, '/x2');
  f = h5read(filename, '/f');
end

assert(lx1==500, 'x1 size')
assert(lx2==1000, 'x2 size')
assert(all(size(f) == [lx1,lx2]), 'array size mismatch')

if ~isinteractive
  return
end
%% PLOT
figure

if (lx2==1)
  plot(x1,f);
  xlabel('x_1')
  ylabel('f')
  title('1-D interp')
else
  imagesc(x2,x1,f);
  axis xy;
  xlabel('x_2')
  ylabel('x_1')
  c=colorbar;
  ylabel(c,'f')
  title('2-D interp')
end
%print -dpng -r300 ~/testinterp.png;

end % function
