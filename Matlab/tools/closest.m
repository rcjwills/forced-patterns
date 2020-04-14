function [index] = closest(x,x0,number,option)

% computed INDEX in vector X closest to X0

% Optional inputs: get only the specified NUMBER of indices if there are several, 
% using OPTION = 'first' or 'last as in find 

if nargin > 2
    index = find(abs(x-x0) == min(abs(x-x0)),number,option);
else
    index = find(abs(x-x0) == min(abs(x-x0)),1,'first');
end