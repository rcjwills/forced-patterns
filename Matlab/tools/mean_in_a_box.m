function [field_box] = mean_in_a_box(lon,lat,field,lat_bnds,lon_bnds,weight)

% area weighted average inside lat_bnds and lon_bnds (2x1 doubles of degree lon/lat limits)
% based on global_mean.m

[x,y] = meshgrid(lon,lat);
if lon_bnds(2) > lon_bnds(1)
    mask = (x >= lon_bnds(1) & x <= lon_bnds(2) & y >= lat_bnds(1) & y <= lat_bnds(2));
else
    mask = ((x <= lon_bnds(2) | x >= lon_bnds(1)) & y >= lat_bnds(1) & y <= lat_bnds(2));
end
nan_mask = ones(size(mask));
nan_mask(mask == 0) = nan;
nan_mask = nan_mask';

if nargin > 5
    if length(size(weight)) > 2
        field_box = global_mean(lon,lat,field.*weight,nan_mask);
    else
        field_box = global_mean(lon,lat,field,nan_mask.*weight);
    end
else
    field_box = global_mean(lon,lat,field,nan_mask);
end