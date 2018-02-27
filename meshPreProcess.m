function [ dynamicMesh ] = meshPreProcess(mesh,porosity,permTensor,NTG)
%TRANCALC Summary of this function goes here
%   Detailed explanation goes here
    if nargin<4
        NTG=ones(1,mesh.nMesh);
    end
    
    dynamicMesh = struct('nMesh',[], 'nFace',[], 'poreVolume0',[], 'depth',[], 'tran',[], 'face2cv',[], 'cv2face',[], 'PVTNUM',[],'SATNUM',[]);
    dynamicMesh.nMesh = mesh.nMesh;
    [nFace,cv2face,face2cv,cv2cv,faceArea,faceNormal,faceCenter]=faceDetection_2D(mesh);
    dynamicMesh.nFace = nFace;
    dynamicMesh.face2cv = face2cv;
    dynamicMesh.cv2face = cv2face;
    dynamicMesh.cv2cv = cv2cv;
    dynamicMesh.cvCenter = mesh.elementCenter;
    dynamicMesh.tran = tranCalc(dynamicMesh,faceArea,faceNormal,faceCenter,permTensor,NTG);
    dynamicMesh.volume = volumeCalc(mesh);
    dynamicMesh.poreVolume0 = poreVolumeCalc(dynamicMesh,porosity,NTG);
    dynamicMesh.depth = depthCalc(mesh);
    dynamicMesh.PVTNUM = ones(1,mesh.nMesh);
    dynamicMesh.SATNUM = ones(1,mesh.nMesh);
    dynamicMesh.maxConnectivity = 4;
end

function tran=tranCalc(mesh,faceArea,faceNormal,faceCenter,permTensor,NTG)
    tran = NaN(1,mesh.nFace);
    CDARCY = 1.0;
    for i=1:mesh.nFace;
        faceCvs = mesh.face2cv(i,:);
        elementCenters = mesh.cvCenter(faceCvs,:);
        gradientVector =  [elementCenters(1,1)-elementCenters(2,1),elementCenters(1,2)-elementCenters(2,2),elementCenters(1,3)-elementCenters(2,3)];
        NTG1 =NTG(faceCvs(1))*[1 1 0]*gradientVector'./norm(gradientVector);
        NTG2 =NTG(faceCvs(2))*[1 1 0]*gradientVector'./norm(gradientVector);
        perm1=NTG1.*permTensor(:,:,faceCvs(1));
        perm2=NTG2.*permTensor(:,:,faceCvs(2));
        perm = 2./(1./perm1+1./perm2);
        normal = faceNormal(i,:);
        tran(i) = CDARCY*(faceArea(i)*abs((perm*gradientVector')'*normal')/(gradientVector*gradientVector'));
    end
end

function poreVolume=poreVolumeCalc(mesh, porosity, NTG)
    poreVolume = mesh.volume .* porosity .*NTG;
end

function depth = depthCalc(mesh)
    cvCoordinate = mesh.elementCenter;
    depth = cvCoordinate(:,3);
end

function [nFace,cv2face,face2cv,cv2cv,faceArea,faceNormal,faceCenter]=faceDetection_2D(mesh)
   
    %face detection for 2D mesh
    nFace=0; cv2cv = NaN(mesh.nMesh,4);
    faceNodes = NaN(mesh.nMesh^2,2);face2cv = NaN(mesh.nMesh^2,2);
    cv2face = NaN(mesh.nMesh,4);nBoundaryFace=0;
    boundaryFaceNodes = NaN(mesh.nMesh^2,2);
    for i=1:mesh.nMesh
        elementNodes = mesh.element2node(i,:);
        for j=1:length(elementNodes)
            nodeElements1 = mesh.node2element(elementNodes(j),:);
            nodeElements1 = nodeElements1(~isnan(nodeElements1));
            if j<length(elementNodes)
                jj=j+1;
            else
                jj=1;
            end
            nodeElements2 = mesh.node2element(elementNodes(jj),:);
            nodeElements2 = nodeElements2(~isnan(nodeElements2));
            internalFace = false;
            for k=1:length(nodeElements1)
                if nodeElements2(nodeElements2 == nodeElements1(k))~= i
                    nFace=nFace+1;
                    faceNodes(nFace,:)=[elementNodes(j), elementNodes(jj)];
                    cv2cv(i,j)= nodeElements1(k);
                    face2cv(nFace,:)=[i,nodeElements1(k)];
                    nans = sum(isnan(cv2face(i,:)));
                    if nans>0, cv2face(i,4-nans+1) = nFace;end
                    internalFace=true;
                    break;
                end
            end
            if ~internalFace
                nBoundaryFace = nBoundaryFace+1;
                boundaryFaceNodes(nBoundaryFace,:)=[elementNodes(j), elementNodes(jj)];
            end
        end
    end
    
    faceNodes = faceNodes(1:nFace,:);face2cv=face2cv(1:nFace,:);
    boundaryFaceNodes = boundaryFaceNodes(1:nBoundaryFace,:);
    faceId =1;
    while faceId <= nFace
        cv1=face2cv(faceId,1);cv2=face2cv(faceId,2);
        i=faceId;
        while i<nFace
            i=i+1;
            cv11=face2cv(i,1);cv22=face2cv(i,2);
            if (cv11 == cv1 && cv22==cv2) || (cv22==cv1 && cv11==cv2)
                face2cv = [face2cv(1:i-1,:);face2cv(i+1:nFace,:)];
                faceNodes = [faceNodes(1:i-1,:);faceNodes(i+1:nFace,:)];
                nFace = nFace-1;
            end
        end
        faceId=faceId+1;
    end
    cvFaceId = ones(mesh.nMesh,1);
    for faceId=1:nFace
        cv1 = face2cv(faceId,1);           cv2 = face2cv(faceId,2);
        cv2face(cv1,cvFaceId(cv1))=faceId; cvFaceId(cv1)=cvFaceId(cv1)+1;
        cv2face(cv2,cvFaceId(cv2))=faceId; cvFaceId(cv2)=cvFaceId(cv2)+1;
    end
    
    faceArea = NaN(1,nFace);faceNormal = NaN(nFace,3);faceCenter=NaN(nFace,3);
    for i=1:nFace
        faceLine=mesh.nodes(faceNodes(i,:),:);
        faceLength = getPolygonLength(faceLine);
        faceArea(i)=faceLength*mesh.increments(3);
        faceNormal(i,:) = [getLinePerpendicular(faceLine,faceLength),0];
        faceCenter(i,:)=[mean(faceLine(:,1)), mean(faceLine(:,2)),0];
    end
end

function volume = volumeCalc(mesh)
    increments = mesh.increments;
    volume = NaN(1,mesh.nMesh);
    for i=1:mesh.nMesh
        volume(i)=getPolygonArea(mesh.nodes(mesh.element2node(i,:),:))*increments(3);
    end
end

