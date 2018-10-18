function validateattr(varargin)
% uses builtin validateattributes() to allow Octave < 4.0 to seamlessly skip validation.
% overloading doesn't work for validateattributes in Octave since it is a core *library* function and there doesn't appear to be a solution besides renaming this function.
v = ver('octave');

if ~isempty(v) && str2double(v.Version) < 4
  return
end

validateattributes(varargin{:})

end