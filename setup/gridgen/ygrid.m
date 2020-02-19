function y = ygrid(ydist, lyp)
narginchk(2,2)
validateattributes(ydist, {'numeric'}, {'scalar','positive'}, mfilename, 'one-way y-distance from origin (meters)',1)
validateattributes(lyp, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'number of y-points',2)

ymin = -ydist/2;
ymax = ydist/2;

if lyp == 1  % degenerate dimension
  % add 2 ghost cells on each side
  y = linspace(ymin, ymax, lyp + 4);
else
  % exclude the ghost cells when setting extents
  y = linspace(ymin, ymax, lyp);
  dy1 = y(2) - y(1);
  dyn = y(end) - y(end-1);
  % now tack on ghost cells so they are outside user-specified region
  y=[y(1)-2*dy1, y(1)-dy1, y, y(end)+dyn, y(end)+2*dyn];
end

end % function