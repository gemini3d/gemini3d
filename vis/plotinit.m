function h = plotinit(xg, visible)
narginchk(1,2)

validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid',1)

if nargin<2,  visible='on'; end
validateattributes(visible, {'char'}, {'vector'}, mfilename, 'figure visibility: on/off', 2)

%Csp = ceil(sqrt(Nt));
%Rsp = ceil(Nt/Csp);

pos2d = [0.1 0.1 0.3 0.3];
pos3d = [0.1 0.1 0.8 0.3];

if any(xg.lx(2:3)==1)  %2D simulation
  figpos = pos2d;
else                            %3D simulation
  figpos = pos3d;
end

if isoctave && strcmp(visible, 'on')
  % octave can have a race condition when making lots of figures leading to
  % "fast fail exception" so we put an arbitrary pause here after EACH
  % figure creation
  % 1 - no crash
  % making arbitrarily large did not help visible plots for Windows.
  % try just saving plots with Octave.
  ts = 1;
else
  ts = 0;
end


h.f1=figure(1);
set(h.f1, 'name', 'V1', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f2=figure(2);
set(h.f2, 'name', 'Ti', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f3=figure(3);
set(h.f3, 'name', 'Te', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f4=figure(4);
set(h.f4, 'name', 'J1', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f5=figure(5);
set(h.f5, 'name', 'V2', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f6=figure(6);
set(h.f6, 'name', 'V3', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f7=figure(7);
set(h.f7, 'name', 'J2', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f8=figure(8);
set(h.f8, 'name', 'J3', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

h.f9=figure(9);
set(h.f9, 'name', 'Topside Potential Phi', 'units', 'normalized', 'position', pos2d, 'visible', visible)
pause(ts)

h.f10=figure(10);
set(h.f10, 'name', 'Ne', 'units', 'normalized', 'position', figpos, 'visible', visible)
pause(ts)

end % function
