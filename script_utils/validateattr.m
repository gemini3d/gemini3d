%!error validateattr([1,1],{'numeric'},{'scalar','positive'})
%!test validateattr(1,{'numeric'},{'scalar','positive'})

function validateattr(varargin)
% uses builtin validateattributes() to allow Octave < 4.0 to seamlessly skip validation.
% overloading doesn't work for validateattributes in Octave since it is a core *library* function % there doesn't appear to be a solution besides renaming this function.
v = ver('octave');

% NOTE: *NOT* isstruct(v) or isfield(v,'Version')
if ~isempty(v) && str2double(v.Version(1:3)) < 4
  return
end

validateattributes(varargin{:})

end
