function [Q] = precip_gaussian2d(pg)
%% makes a 2D Gaussian shape in Latitude, Longitude
%
narginchk(1, 1)
validateattributes(pg, {'struct'}, {'scalar'}, mfilename, "precip grid", 1)

Q = 10*exp(-(pg.MLON - pg.mlon_mean).^2 / (2*pg.mlon_sigma^2)) ...
    .* exp(-(pg.MLAT - pg.mlat_mean).^2 / (2*pg.mlat_sigma^2));

end % function