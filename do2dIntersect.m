function [ output ] = do2dIntersect( line1, line2 )
%DO2DINTERSECT Summary of this function goes here
%   Detailed explanation goes here
    p1 = line1(1,:);
    q1 = line1(2,:);
    p2 = line2(1,:);
    q2 = line2(2,:);
    
    orient1 = get2dOrientation(p1,q1,p2);
    orient2 = get2dOrientation(p1,q1,q2);
    orient3 = get2dOrientation(p2,q2,p1);
    orient4 = get2dOrientation(p2,q2,q1);
    
    output_General  = (orient1 ~= orient2) && (orient3 ~= orient4);
    output_Exception= (orient1 == 0)&&(isOnSegment(p2,[p1;q1]));
    output_Exception= output_Exception || (orient2 == 0)&&(isOnSegment(q2,[p1;q1]));
    output_Exception= output_Exception || (orient3 == 0)&&(isOnSegment(p1,[p2;q2]));
    output_Exception= output_Exception || (orient4 == 0)&&(isOnSegment(q1,[p2;q2]));
    
    
    output = output_General || output_Exception;

end

