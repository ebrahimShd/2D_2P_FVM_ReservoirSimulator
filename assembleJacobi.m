function [ jacobi ] = assembleJacobi( mesh )
    nPhases = mesh.nPhases;
    nCv = mesh.nMesh;
    jacobi = zeros(nPhases*nCv,nPhases*nCv);
    nUpwind = mesh.nUpwind;
    
    for cvId = 1:nCv
        jacobi(nPhases*cvId-1:nPhases*cvId,nPhases*cvId-1:nPhases*cvId)=[mesh.temporalDerivatives(cvId,1:2);mesh.temporalDerivatives(cvId,3:4)]-[mesh.wellDerivatives(cvId,1:2);mesh.wellDerivatives(cvId,3:4)];
    end
    
    isUp = zeros(nUpwind,mesh.nFace,2);
    for i=1:nUpwind
        isUp(i,:,:) = (mesh.face2cv==(mesh.upWind(:,i)*[1 1]));
    end

    for faceId = 1:mesh.nFace
        faceTran = mesh.tran(faceId);
        faceCvs = mesh.face2cv(faceId,:);
        cvUps = mesh.upWind(faceId,:);
        for eq = 1:nPhases
            mobility = mesh.mobility(faceId,eq);
            dp = mesh.deltaP(faceId,eq);
            for var = 1:nPhases
                dTrd = mesh.dTrd{2*(eq-1)+var};     
                dDpd = mesh.dDpd{2*(eq-1)+var};
                for faceSide=1:2
                    cvId = faceCvs(faceSide);
                    row = nPhases*(cvId-1)+eq;
                    deltaP = (-2*(faceSide-1)+1)*dp;
                    for toCv=1:2
                        toCvId = faceCvs(toCv);
                        column = nPhases*(toCvId-1)+var;
                        dDeltaP = (-2*(faceSide-1)+1)*dDpd(faceId,toCv);
%                         [isUp,index] = isinUpwinds(mesh,toCvId,faceId);
%                         if isUp,facedTrd = dTrd(faceId,index);else facedTrd = 0.0;end
                        index = find(isUp(:,faceId,toCv),1);
                        if ~isempty(index),facedTrd = dTrd(faceId,index);else facedTrd = 0.0;end
                        jacobi(row,column) = jacobi(row,column)-faceTran*(facedTrd*deltaP+mobility*dDeltaP);
                    end
                    if nUpwind>1
                        for index = 1:nUpwind
                            if isempty(find(isUp(index,faceId,:),1))
                                upId = cvUps(index);
                                dDeltaP = 0.0;
                                facedTrd = dTrd(faceId,index);
                                column = nPhases*(upId-1)+var;
                                jacobi(row,column) = jacobi(row,column)-faceTran*(facedTrd*deltaP+mobility*dDeltaP);
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    jacobi = sparse(jacobi);
end

function [out,index]=isinUpwinds(mesh,cvId, faceId)
    cvUps = mesh.upWind(faceId,:);
    index = find(cvUps == cvId);
    out = ~isempty(index);
end

function out = isInFaceSide(mesh,upId,faceId)
    out = false;
    faceCvs = mesh.face2cv(faceId,:);
    %out = ~isempty(find(faceCvs==upId,1));
    for i=1:2
        if upId==faceCvs(i)
            out=1;
            break;
        end
    end
end