function [ val ] = findSValue( X,Y,letter,S )
%FINDSVALUE Summary of this function goes here
%   Detailed explanation goes here
val=S(letter==X,letter==Y);

end

