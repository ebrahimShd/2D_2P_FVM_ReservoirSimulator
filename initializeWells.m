function [ wells ] = initializeWells( numberOfWells, wellPositions, geoMesh, dynaMesh, wellRws, types, phases, constarints, constraintValues,wellIndices)
    wells = struct('nWells',[], 'well2cv', [], 'wellIndex', [],'cv2well',[], 'wellType',[],'wellPhase',[],'activeConstraint', [], 'constraintValue', []);
    well2cv=NaN(1,numberOfWells);cv2well=NaN(1,geoMesh.nMesh);
    for wellId=1:numberOfWells
        for cvId=1:geoMesh.nMesh
            cvNodes = geoMesh.nodes(geoMesh.element2node(cvId,:),:);
            cvNodes = [cvNodes;cvNodes(1,:)];
            if isInside(wellPositions(wellId,:),cvNodes(:,1:2));
                well2cv(wellId)=cvId; 
                cv2well(cvId)=wellId; break;
            end
        end
    end
    
    
    if nargin<5, wellIndices = getWI(numberOfWells, well2cv, geoMesh, dynaMesh, wellRws);end
    
    wells.nWells = numberOfWells;
    wells.well2cv = well2cv;
    wells.wellIndex = wellIndices;
    wells.cv2well = cv2well;
    wells.wellType = types;
    wells.wellPhase = phases;
    wells.activeConstraint = constarints;
    wells.constraintValue = constraintValues;

end

function wellIndices=getWI(nWells, well2cv, geoMesh, dynaMesh, wellRws)
    wellIndices = 1*ones(1,numberOfWells); 
        for wellId=1:nWells
            cvId = well2cv(wellId);
            dx = geoMesh.increments(cvId,1);
            dy = geoMesh.increments(cvId,2);
            dz = geoMesh.increments(cvId,3);
            
            
        end
end
