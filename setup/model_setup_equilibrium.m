function model_setup_equilibrium(p)
%% setup equilibrium simulation
% this is to be called by model_setup.m

arguments
  p (1,1) struct
end

%% GRID GENERATION
xg = makegrid_cart_3D(p);

writegrid(p, xg);

%% Equilibrium input generation

[ns,Ts,vsx1] = eqICs3D(p, xg);
% Note: should be rewritten to include the neutral module form the fortran code
writedata(p.ymd, p.UTsec0, ns, vsx1, Ts, p.simdir, p.format, p.realbits);

end % function
