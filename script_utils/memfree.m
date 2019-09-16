%!assert(memfree() > 0)
%!assert(isnumeric(memfree))
function freebytes = memfree()
%% find free physical RAM on Windows (with or without Cygwin) and Linux systems
% currently Matlab doesn't support memory() on Linux/Mac systems
% This function is meant to give free memory using Matlab or Octave
%
% It demonstrates using Python from Matlab or GNU Octave seamlessly.
%
% Output:
% --------
% free physical RAM [bytes]
%
% If Python psutils not available, returns NaN
%
% Michael Hirsch, Ph.D.

try
  freebytes = double(py.psutil.virtual_memory().available);
catch
  [~,freebytes] = system('python -c "import psutil; print(psutil.virtual_memory().available)"');
  freebytes = str2double(freebytes);
end

end %function
