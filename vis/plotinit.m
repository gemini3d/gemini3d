function h = plotinit(xg)

if nargin<1, xg=[]; end

%Csp = ceil(sqrt(Nt));
%Rsp = ceil(Nt/Csp);

h.f1=figure('name','V1', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f2=figure('name','Ti', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f3=figure('name','Te', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f4=figure('name','J1', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f5=figure('name','V2', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f6=figure('name','V3', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f7=figure('name','J2', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
h.f8=figure('name','J3', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
if ~isempty(xg) && xg.lx(2)>1 && xg.lx(3)>1 % a 3-D simulation
  h.f9=figure('name','phiTop', 'units', 'normalized', 'position', [.1, .1, .5, .5]);
else
  h.f9 = [];
end
h.f10=figure('name','Ne', 'units', 'normalized', 'position', [.1, .1, .5, .5]);

end
