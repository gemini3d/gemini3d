function close_enough = allclose(actual, desired, rtol, atol)
% close_enough = allclose(actual, desired, rtol, atol)
%
% Inputs
% ------
% atol: absolute tolerance (scalar)
% rtol: relative tolerance (scalar)
%
% Output
% ------
% close_enough: logical TRUE if "actual" is close enough to "desired"
%
% based on numpy.testing.assert_allclose
% https://github.com/numpy/numpy/blob/v1.13.0/numpy/core/numeric.py#L2522
% for Matlab and GNU Octave
%
% if "actual" is within atol OR rtol of "desired", return true

narginchk(2,4)
validateattr(actual, {'numeric'}, {'nonempty'}, mfilename, 'measured values', 1)
validateattr(desired, {'numeric'}, {'nonempty'}, mfilename, 'desired reference values', 2)
if nargin < 3 || isempty(rtol)
  rtol=1e-8;
else
  validateattr(rtol, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'relative tolerance', 3)
end
if nargin < 4 || isempty(atol)
  atol = 1e-9;
else
  validateattr(atol, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'absolute tolerance', 4)
end

%% compare
actual = actual(:);
desired = desired(:);

measdiff = abs(actual-desired);
tol = atol + rtol * abs(desired);

close_enough = all(measdiff <= tol);

end  % function