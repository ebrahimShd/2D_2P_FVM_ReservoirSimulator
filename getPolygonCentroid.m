function [ centroid ] = getPolygonCentroid( polygon, area )
%GETPOLYGONCENTROID Summary of this function goes here
%   Detailed explanation goes here
    if nargin<2
        area = getPolygonArea(polygon);
    end
    x = 0;y=0;xs=polygon(:,1);ys=polygon(:,2);
    for i=1:size(polygon,1)-1
        x = x + (xs(i) + xs(i+1))*(xs(i)*ys(i+1)-xs(i+1)*ys(i));
        y = y + (ys(i) + ys(i+1))*(xs(i)*ys(i+1)-xs(i+1)*ys(i));
    end
    centroid = 1/(6*area)*[x,y];
end

