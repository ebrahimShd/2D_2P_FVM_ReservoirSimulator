function [ mesh ] = meshGenerator2D( boundary, meshSize, dz, dx, dy)

    %dimention of problem
    nI = meshSize(1); nJ = meshSize(2); nMesh = nI*nJ; nNodes = (nI+1)*(nJ+1);
    
    %variable initialization
    mesh = struct('nMesh',[], 'nNode',[], 'nodes',[], 'elementCenter',[],'increments',[], 'node2element',[], 'element2node',[]);
    
    nodes = NaN(nNodes,3); centerPoints = NaN(nMesh,3);
    element2node = NaN(nMesh,4); node2element = NaN(nNodes,4);
    
    isElementActive = true(1,nMesh); isNodeActive = true(1,nNodes);
    nElementActive = nMesh; nNodeActive = nNodes;
    
    %geometry calculation
    xMin = min(boundary(:,1)); xMax = max(boundary(:,1));
    yMin = min(boundary(:,2)); yMax = max(boundary(:,2));
    if nargin<4
        dx=(xMax-xMin)/nI*ones(1,nI); 
        dy=(yMax-yMin)/nJ*ones(1,nJ);
    elseif nargin<5
        dy=(yMax-yMin)/nJ*ones(1,nJ);
    end
    if abs(sum(dx)-(xMax-xMin))>1e-5,dx = (dx./sum(dx)).*(xMax-xMin);end
    if abs(sum(dy)-(yMax-yMin))>1e-5,dy = (dy./sum(dy)).*(yMax-yMin);end
    xs = NaN(1,nI+1); ys = NaN(1,nJ+1);xs(1)=xMin;ys(1)=yMin;
    for i=2:nI+1,xs(i)=xs(i-1)+dx(i-1);end
    for i=2:nJ+1,ys(i)=ys(i-1)+dy(i-1);end
    
    %nodes to element mapping
    id = 1;
    for j=1:nJ+1
        for i=1:nI+1
            nodes(id,:)=[xs(i),ys(j),0];
            element1=NaN;element2=NaN;element3=NaN;element4=NaN;
            if i>1 && j<nJ+1, element3=(j-1)*nI+i-1;end
            if j>1 && i<nI+1, element2=(j-2)*nI+i;end
            if i>1 && j>1, element1=(j-2)*nI+i-1;end
            if i<nI+1 && j<nJ+1, element4=(j-1)*nI+i;end
            node2element(id,:)=[element1, element2, element3, element4];
            id=id+1;
        end
    end
    
    %detect inactive elements
    id = 1;
    for j=1:nJ
        for i=1:nI
            element2node(id,:)=[(j-1)*(nI+1)+i, (j-1)*(nI+1)+i+1, j*(nI+1)+i+1, j*(nI+1)+i];
            elementNodes = nodes(element2node(id,:),:);
            centerPoints(id,:)=[mean(elementNodes(:,1)),mean(elementNodes(:,2)),0];
            if ~isInside(centerPoints(id,1:2),boundary)
                isElementActive(id)=false;
                nElementActive = nElementActive - 1;
            end
            id = id + 1;
        end
    end
    %detect inactive nodes
    for id=1:nNodes
        elements = node2element(id,:); elements = elements(~isnan(elements));
        deActivate = true;
        for j=1:length(elements)
            if isElementActive(elements(j))
                deActivate = false;
            end
        end
        if deActivate
            isNodeActive(id)=false;
            nNodeActive = nNodeActive-1;
        end
    end
    
    % element inactive 2 active mapping
    id=1;element2active = NaN(1,nMesh);active2element=NaN(1,nElementActive);
    for i=1:nMesh
        if isElementActive(i)
            element2active(i)=id;
            active2element(id)=i;
            id = id+1;
        end
    end
    
    % element inactive 2 active mapping
    id=1;node2active=NaN(1,nNodes);active2node=NaN(1,nNodeActive);
    for i=1:nNodes
        if isNodeActive(i)
            node2active(i)=id;
            active2node(id)=i;
            id = id+1;
        end
    end
    
    % evaluate element prop for active only elements
    activeNode2element=NaN(nNodeActive,4);
    activeCenterPoints = centerPoints(active2element,:);
    activeElement2node = node2active(element2node(active2element,:));
    
    for i=1:nNodeActive
        elements = node2element(active2node(i),:);elements = elements(~isnan(elements));
        elements = element2active(elements);elements = [elements,NaN(1,4-length(elements))];
        activeNode2element(i,:)=elements;
    end
    
    mesh.nMesh = nElementActive;
    mesh.nNode = nNodeActive;
    mesh.nodes = nodes(active2node,:);
    mesh.increments = [mean(dx),mean(dy),dz];
    mesh.elementCenter = activeCenterPoints;
    mesh.node2element = activeNode2element;
    mesh.element2node = activeElement2node;
end

