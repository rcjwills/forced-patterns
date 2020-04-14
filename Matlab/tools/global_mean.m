function [field_global_mean,field_departure] = global_mean(lon,lat,field,weight)

% global mean of field on lat/lon grid defined in input lat and lon
% use weight to density weight for non vertically integrated fields
% assumes lon & lat are the first 2 dimensions

lat = lat.*pi./180;
dlat = zeros(1,length(lat));
for j = 1:length(lat)
    if j == 1
        dlat(j) = abs(lat(2)-lat(1));
    elseif j == length(lat)
        dlat(j) = abs(lat(j)-lat(j-1));
    else
        dlat(j) = abs((lat(j+1)-lat(j-1)))/2;
    end
end

    dlat = repmat(dlat,length(lon),1);
    s = size(lat);
    if s(1) > 1
        lat = lat';
    end
    lat = repmat(lat,length(lon),1);
    dlon = abs(lon(2)-lon(1))*pi/180;
    
    s = size(field);
    
    ones_nans = ones(size(field));
    ones_nans(isnan(field)) = nan;
    
    if length(s) > 2 % 3D or 4D input
        if length(s) > 3 % 4D input (e.g. lon, lat, sig, time)
            dlat = repmat(dlat,[1 1 s(3) s(4)]);
            lat = repmat(lat,[1 1 s(3) s(4)]);
            if nargin > 3
                weight = repmat(weight,[1 1 s(3) s(4)]);
                weighted_field = reshape(dlat.*dlon.*cos(lat).*field.*weight,[s(1)*s(2),s(3),s(4)]);
                weight_vector = reshape(dlat.*dlon.*cos(lat).*ones_nans.*weight,[s(1)*s(2),s(3),s(4)]);
            else
                weighted_field = reshape(dlat.*dlon.*cos(lat).*field,[s(1)*s(2),s(3),s(4)]);
                weight_vector = reshape(dlat.*dlon.*cos(lat).*ones_nans,[s(1)*s(2),s(3),s(4)]);
            end
        else
            dlat = repmat(dlat,[1 1 s(3)]);
            lat = repmat(lat,[1 1 s(3)]);
            if nargin > 3 % weighted global mean
                weight = repmat(weight,[1 1 s(3)]);
                weighted_field = reshape(dlat.*dlon.*cos(lat).*field.*weight,[s(1)*s(2),s(3)]);
                weight_vector = reshape(dlat.*dlon.*cos(lat).*ones_nans.*weight,[s(1)*s(2),s(3)]);
            else
                weighted_field = reshape(dlat.*dlon.*cos(lat).*field,[s(1)*s(2),s(3)]);
                weight_vector = reshape(dlat.*dlon.*cos(lat).*ones_nans,[s(1)*s(2),s(3)]);
            end
        end
    else % 2D input (lon, lat)
        if nargin > 3 % weighted global mean
            weighted_field = reshape(dlat.*dlon.*cos(lat).*field.*weight,[s(1)*s(2),1]);
            weight_vector = reshape(dlat.*dlon.*cos(lat).*ones_nans.*weight,[s(1)*s(2),1]);
        else
            weighted_field = reshape(dlat.*dlon.*cos(lat).*field,[s(1)*s(2),1]);
            weight_vector = reshape(dlat.*dlon.*cos(lat).*ones_nans,[s(1)*s(2),1]);
        end
    end
 
    field_global_mean = nanmean(weighted_field,1);
    ones_global_mean = nanmean(weight_vector,1);
    field_global_mean = field_global_mean./ones_global_mean;

try
if size(field,4) > 1
    field_departure = field - repmat(reshape(field_global_mean,[1 1 s(3) s(4)]),[s(1) s(2) 1 1]);
else
    if size(field_global_mean,3) > 1
        field_departure = field - repmat(field_global_mean,[s(1) s(2) 1]);
    else
        field_departure = field - field_global_mean;
    end
end
catch
    if nargout > 1
        disp('error with departure from global mean calculation')
        field_departure = nan;
    end
end

field_global_mean = squeeze(field_global_mean);

