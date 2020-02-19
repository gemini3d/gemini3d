function x = xgrid(xdist, lxp)
narginchk(2,2)
validateattributes(xdist, {'numeric'}, {'scalar','positive'}, mfilename, 'one-way x-distance from origin (meters)',1)
validateattributes(lxp, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'number of x-points',2)

xmin = -xdist/2;
xmax = xdist/2;

if lxp == 1  % degenerate dimension
  % add 2 ghost cells on each side
  x = linspace(xmin, xmax, lxp + 4);
else
  % exclude the ghost cells when setting extents
  x = linspace(xmin, xmax, lxp);
  dx1 = x(2) - x(1);
  dxn = x(end) - x(end-1);
  % now tack on ghost cells so they are outside user-specified region
  x=[x(1)-2*dx1, x(1)-dx1, x, x(end)+dxn, x(end)+2*dxn];
end

end % function