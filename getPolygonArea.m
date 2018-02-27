function [ area ] = getPolygonArea( polygon )
%GETPOLYGONAREA Summary of this function goes here
%   Detailed explanation goes here
    area=0;
    nPoints = size(polygon,1);
    if ~isequal(polygon(1,:),polygon(nPoints,:))
        polygon(nPoints+1,:)=polygon(1,:);
        nPoints = nPoints+1;
    end
    for i=1:nPoints-1
        area = area + polygon(i,1)*polygon(i+1,2)-polygon(i+1,1)*polygon(i,2);
    end
    area = area/2;
end

