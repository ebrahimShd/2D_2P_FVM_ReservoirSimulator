function [ output ] = isOnSegment( point, segment )
%ISONSEGMENT Summary of this function goes here
%   Detailed explanation goes here
    xp = point(1);
    yp = point(2);
    xSegment = segment(:,1)';
    ySegment = segment(:,2)';
    output = (xp<=max(xSegment) || xp>=min(xSegment))&&((yp<=max(ySegment))&&(yp>=min(ySegment)));
end

