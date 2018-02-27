function [ output ] = isInside( point, polygon)
%ISINSIDE Summary of this function goes here
%   Detailed explanation goes here
    INF = max(polygon(:,1))*10;
    exLine = [point;INF point(2)];
    count = 0;
    onVertex = false;
    for i=1:size(polygon,1)-1
        count = count + do2dIntersect(exLine, [polygon(i,:);polygon(i+1,:)]);
        onVertex = onVertex || isOnSegment(point,[polygon(i,:),polygon(i+1,:)]);
    end
    for i=1:size(polygon,1)
        if isequal(point,polygon(i))
            onVertex = true;
            break;
        end
    end
    output = mod(count,2)>0;
    output = onVertex || output;
    

end

