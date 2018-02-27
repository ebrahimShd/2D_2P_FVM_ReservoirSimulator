function [ vector ] = getLinePerpendicular( line, length )
%GETVECTORPERPENDICULAR Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 2
        length = getPolygonLength(line);
    end
    vector = [line(1,1)-line(2,1),line(1,2)-line(2,2)];
    vector = [vector(2),-vector(1)] ./ length;
end

