function pg = precip_grid(xg, p, llat, llon)

narginchk(4, 4)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, "spatial grid", 1)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, "params", 2)
validateattributes(llat, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'number of latitude', 3)
validateattributes(llat, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'number of longitude', 4)

thetamin = min(xg.theta(:));
thetamax = max(xg.theta(:));
mlatmin = 90-thetamax*180/pi;
mlatmax = 90-thetamin*180/pi;
mlonmin = min(xg.phi(:))*180/pi;
mlonmax = max(xg.phi(:))*180/pi;

% add a 1% buffer
latbuf = 1/100*(mlatmax-mlatmin);
lonbuf = 1/100*(mlonmax-mlonmin);
% create the lat,lon grid
pg.mlat = linspace(mlatmin-latbuf,mlatmax+latbuf,llat);
pg.mlon = linspace(mlonmin-lonbuf,mlonmax+lonbuf,llon);
[pg.MLON, pg.MLAT] = ndgrid(pg.mlon, pg.mlat);
pg.mlon_mean = mean(pg.mlon);
pg.mlat_mean = mean(pg.mlat);

%% disturbance width
mlatsig = p.precip_latwidth*(mlatmax-mlatmin);
% to avoid divide by zero below
pg.mlat_sigma = max(mlatsig, 0.01);
pg.mlon_sigma = p.precip_lonwidth*(mlonmax-mlonmin);

end % function
