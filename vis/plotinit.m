function h = plotinit(xg, visible)

if nargin<1
  xg=[]; 
elseif ~isempty(xg)
  validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid',1)
end

if nargin<2
  visible='on';
else
  validateattributes(visible, {'char'}, {'vector'}, mfilename, 'figure visibility', 2)
end

%Csp = ceil(sqrt(Nt));
%Rsp = ceil(Nt/Csp);

%figpos=[0.1 0.1 0.5 0.5];
figpos=[0.1 0.1 0.8 0.3];

h.f1=figure('name','V1', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f2=figure('name','Ti', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f3=figure('name','Te', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f4=figure('name','J1', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f5=figure('name','V2', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f6=figure('name','V3', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f7=figure('name','J2', 'units', 'normalized', 'position', figpos, 'visible', visible);
h.f8=figure('name','J3', 'units', 'normalized', 'position', figpos, 'visible', visible);
if ~isempty(xg) && xg.lx(2)>1 && xg.lx(3)>1 % a 3-D simulation
  h.f9=figure('name','phiTop', 'units', 'normalized', 'position', figpos, 'visible', visible);
else
  h.f9 = [];
end
h.f10=figure('name','Ne', 'units', 'normalized', 'position', figpos, 'visible', visible);

end
