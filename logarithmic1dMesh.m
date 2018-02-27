function [ dxs ] = logarithmic1dMesh( interval, central, nx, dxmin)
%LOGARITHMIC1DMESH Summary of this function goes here
%   Detailed explanation goes here
    if central
        domainDx = interval(2)-interval(1);
        if mod(nx,2)==0,nx=nx+1;end
        nLog = (nx-1)/2;
        xmid=interval(1)+domainDx/2; xmin=xmid+dxmin/2; xmax=interval(2);
        xmax = xmax-xmin; xmin = dxmin;
        xminLog = log(xmin); xmaxLog = log(xmax);
        hLog = (xmaxLog-xminLog)/(nLog-1);
        xLog = xminLog:hLog:xmaxLog;
        xs =  exp(xLog) + (xmid+dxmin/2)*ones(size(xLog)) ;
        xs = [xmid+dxmin/2 xs];
        xleft = zeros(1,nx/2+1/2);
        for i=(nx+1)/2:-1:1
            xleft(i) = xmid - (xs(12-i)-xmid);
        end
        xs = [xleft xs];
    end
    dxs = zeros(1,nx);
    for i=1:nx
        dxs(i)=xs(i+1)-xs(i);
    end
end

