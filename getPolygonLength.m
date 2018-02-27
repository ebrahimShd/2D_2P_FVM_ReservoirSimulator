function [ length ] = getPolygonLength( polygon )
%GETPOLYGONLENGTH Summary of this function goes here
%   Detailed explanation goes here
    length=0;
    for i=1:size(polygon,1)-1
        length = length + sqrt((polygon(i,1)-polygon(i+1,1))^2+(polygon(i,2)-polygon(i+1,2))^2);
    end

end

