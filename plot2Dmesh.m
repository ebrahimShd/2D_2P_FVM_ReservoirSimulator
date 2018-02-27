function [ ] = plot2Dmesh( mesh )
%PLOT2DMESH Summary of this function goes here
%   Detailed explanation goes here
    figure(1);
    hold on;
    for i=1:mesh.nMesh
        cvNodes = mesh.nodes(mesh.element2node(i,:),:);
        for j=1:length(cvNodes)-1
            plot(cvNodes(j:j+1,1),cvNodes(j:j+1,2),'*-b');
        end
    end
    plot(mesh.elementCenter(:,1),mesh.elementCenter(:,2),'or');
end

