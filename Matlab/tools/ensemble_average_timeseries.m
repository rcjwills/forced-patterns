function [Xe] = ensemble_average_timeseries(X,ne)

l = size(X,1)/ne;
Xe = zeros([l size(X,2)]);

for i = 1:ne
    i1 = (i-1)*l+1;
    i2 = i*l;
    Xe = Xe + X(i1:i2,:)./ne;
end
