function h = plotinit(xg, visible)

if nargin<1
  xg=[]; 
elseif ~isempty(xg)
  validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid',1)
end

if nargin<2
  visible='on';
else
  validateattributes(visible, {'char'}, {'vector'}, mfilename, 'figure visibility: on/off', 2)
end

%Csp = ceil(sqrt(Nt));
%Rsp = ceil(Nt/Csp);

if(xg.lx(3)==1 || xg.lx(2)==1)  %2D simulation
  figpos=[0.1 0.1 0.3 0.3];
else                            %3D simulation
  figpos=[0.1 0.1 0.8 0.3];    
end


h.f1=figure(1);
set(h.f1, 'name', 'V1', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f2=figure(2);
set(h.f2, 'name', 'Ti', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f3=figure(3);
set(h.f3, 'name', 'Te', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f4=figure(4);
set(h.f4, 'name', 'J1', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f5=figure(5);
set(h.f5, 'name', 'V2', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f6=figure(6);
set(h.f6, 'name', 'V3', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f7=figure(7);
set(h.f7, 'name', 'J2', 'units', 'normalized', 'position', figpos, 'visible', visible)
h.f8=figure(8);
set(h.f8, 'name', 'J3', 'units', 'normalized', 'position', figpos, 'visible', visible)
if ~isempty(xg) && xg.lx(2)>1 && xg.lx(3)>1 % a 3-D simulation
  h.f9=figure(9);
  set(h.f9, 'name','phiTop', 'units', 'normalized', 'position', figpos, 'visible', visible)
else
  h.f9 = [];
end
h.f10=figure(10);
set(h.f10, 'name', 'Ne', 'units', 'normalized', 'position', figpos, 'visible', visible)

end
