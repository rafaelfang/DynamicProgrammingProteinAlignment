function [ maxVal ] = getMax( varargin )
%find and return the max value from numbers of input
%   written by Chao Fang



    maxVal=-999999;
    nVarargs = length(varargin);

    for k = 1:nVarargs
        maxVal=max(maxVal,varargin{k});
    end
end

