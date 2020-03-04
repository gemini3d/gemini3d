function model_setup_interp(p)
%% setup interpolated simulation based on equilibrium simulation
% this is to be called by model_setup.m

arguments
  p (1,1) struct
end
%% GRID GENERATION
xg = makegrid_cart_3D(p);

eq2dist(p, xg);

%% potential boundary conditions
if isfield(p, 'flagE0file') && p.flagE0file
  if p.lxp == 1 || p.lyp == 1
    Efield_BCs_2d(p, xg);
  else % 3D
    Efield_BCs_3d(p, xg);
  end
end

%% aurora
if isfield(p, 'flagprecfile') && p.flagprecfile
  particles_BCs(p, xg)
end

end % function
