function expanded = expanduser(p)
%%
% For now, handles only a leading tilde, does not currently handle specifying ~otheruser
%
% isunix==1 (linux & cygwin) example:
% expanduser('~/Downloads/foo')
% ans = /home/joespc/Downloads/foo
%
% ispc==1 example (Windows)
% expanduser('~/Downloads/foo')
% ans = C:\joespc/Downloads/foo
%
% Useful for Matlab functions like h5read() and some Computer Vision toolbox functions
% that can't handle ~ and Matlab does not consider it a bug per conversations with 
% Mathworks staff
%
% Michael Hirsch
%
% tested with Matlab and Octave on Windows, Cygwin, Linux, and WINE
%
%% try python first
try %requires Matlab R2014b or newer for following line
    expanded = char(py.os.path.expanduser(p));
    return
end
%% if you have old Matlab or Octave
%% what is the home path
if ispc % windows
    home = [getenv('HOMEDRIVE'),getenv('HOMEPATH')];
else %linux,mac
    home = getenv('HOME');
end %if

if isempty(home)
    warning('empty HOME environment variable, returning unmodified path')
    expanded =p;
    return
end %if
%% now let's look at your path, does it have a leading tilde?
if ~isempty(p) && ischar(p) && size(p,1) == 1
    if length(p) == 1 && strcmp(p,'~')
        expanded = home;
    elseif strcmp(p(1:2),'~/') || strcmp(p(1:2),'~\')
        expanded = [home,p(2:end)];
    elseif ~isempty(regexp(p,'~.*/', 'once')) || ~isempty(regexp(p,'~.*\\', 'once'))
            warning('the ~otheruser case is not handled yet')
            expanded = p;
    else
        expanded = p; 
    end %if
else
    warning('i only handle non-array strings for now') %TODO: consider cellfun()
    expanded = p;
end %if

end %function
