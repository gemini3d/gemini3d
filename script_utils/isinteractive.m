%!assert(islogical(isinteractive))
function isinter = isinteractive()
 %% tell if the program is being run interactively or not.
 % helpful to say pause after making groups of plots--only if user has GUI desktop open.
 % don't use batchStartupOptionUsed as it neglects the "-nodesktop" case
persistent inter;

if isempty(inter)
    if isoctave
      inter = isguirunning;
    else % matlab, this test below doesn't work for Octave
	  inter = usejava('desktop');
    end
end

isinter=inter; % has to be a separate line/variable for matlab

end
