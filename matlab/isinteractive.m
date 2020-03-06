%!assert(islogical(isinteractive))
function isinter = isinteractive()
%% tell if the program is being run interactively or not.

persistent inter;

if isempty(inter)
  if isoctave
    inter = isguirunning;
  else
    % matlab, this test doesn't work for Octave
    % don't use batchStartupOptionUsed as it neglects the "-nodesktop" case
    inter = usejava('desktop');
  end
end

% has to be a separate line / variable for matlab
isinter=inter;

end
