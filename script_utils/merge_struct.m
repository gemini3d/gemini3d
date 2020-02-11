%!assert(isstruct(merge_struct(struct(), struct())))
function s1 = merge_struct(s1, s2, overwrite)
%% merge_struct(struct, s)
%
% merges fields of SCALAR structs "s1" and "s2", optionally overwriting existing
% s1 fields from s2.

narginchk(2, 3)

validateattributes(s1, {'struct'}, {'scalar'})
validateattributes(s2, {'struct'}, {'scalar'})
if nargin < 3
  overwrite = false;
end
validateattributes(overwrite, {'logical'}, {'scalar'})

s1fields = fieldnames(s1);
for field = fieldnames(s2)'
  if ~overwrite && isfield(s1fields, field{:})
    error(['not overwriting field ', field{:}])
  end
  s1.(field{:}) = s2.(field{:});
end
end % function