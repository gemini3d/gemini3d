function [flag, suffix] = printflag(fmt)
narginchk(1,1)
validateattributes(fmt, {'char'}, {'vector'}, mfilename, 'png or eps', 1)

flag = []; suffix = [];

switch fmt
  case 'png', flag = '-dpng'; suffix = '.png';
  case 'eps', flag = '-depsc2'; suffix = '.eps';
end

end
