function residuals=evaluateResiduals(mesh,well,vectorized)
    global msh wel
    msh = mesh;wel = well;
    nCv = msh.nMesh;
    nPhases = msh.nPhases;
    if vectorized
        temp = vectorizedTemporal();
        spac = vectorizedSpacial();
        well = vectorizedWell();
        residuals = temp-spac-well;
        residuals = reshape(residuals',[nCv*nPhases,1]);
    else
        nPhases = mesh.nPhases;
        residuals = zeros(nPhases*mesh.nMesh,1);
        for cvId = 1:mesh.nMesh
            temp = evaluateTemporal(cvId);
            spac = evaluateSpacial(cvId);
            well = evaluateWell(cvId);
            for phase = 1:nPhases
                residuals(nPhases*(cvId-1)+phase) = temp(phase)-spac(phase)-well(phase);
            end
        end
    end
end

function output = evaluateTemporal(cvId)
    global msh
    output = zeros(1,2);
    sw = msh.sw(cvId);      swPrev = msh.swPrev(cvId);
    so = msh.so(cvId);      soPrev = msh.soPrev(cvId);
    bo = msh.bo(cvId);      boPrev = msh.boPrev(cvId);
    bw = msh.bw(cvId);      bwPrev = msh.bwPrev(cvId);
    pv = msh.poreVolume(cvId);    pvPrev = msh.pvPrev(cvId);
    dt = msh.dt;
    output(1) = (pv*so/bo-pvPrev*soPrev/boPrev)/dt;
    output(2) = (pv*sw/bw-pvPrev*swPrev/bwPrev)/dt;
end

function output = vectorizedTemporal()
    global msh
    dt = msh.dt;
    output = [(msh.poreVolume.*msh.so./msh.bo-msh.pvPrev.*msh.soPrev./msh.boPrev)./dt;(msh.poreVolume.*msh.sw./msh.bw-msh.pvPrev.*msh.swPrev./msh.bwPrev)./dt]';
end

function output = evaluateSpacial(cvId)
    global msh
    output = zeros(1,2);
    cvFaces = msh.cv2face(cvId,:);cvFaces = cvFaces(~isnan(cvFaces)); nFace = length(cvFaces);
    for i = 1:nFace
		faceId = cvFaces(i);
        tran = msh.tran(faceId);
        faceCvs = msh.face2cv(faceId,:);
        mobility = msh.mobility(faceId,:);
        dp = msh.deltaP(faceId,:); if cvId == faceCvs(2), dp = -1.*dp;end
        for ii=1:2
            output(ii) = output(ii)+tran*mobility(ii)*dp(ii);
        end
    end
end

function output = vectorizedSpacial()
    global msh
    outputFace=msh.tran'*ones(1,2).*msh.mobility.*msh.deltaP;
    output = zeros(msh.nMesh,2);
    for cvId = 1:msh.nMesh
        cvFaces = msh.cv2face(cvId,:);cvFaces = cvFaces(~isnan(cvFaces));
        faceCvs = msh.face2cv(cvFaces,:);
        dpMult = 2*(faceCvs(:,1)==cvId)-ones(size(faceCvs,1),1);
        output(cvId,:) = sum(dpMult*ones(1,2).*outputFace(cvFaces,:),1);
    end
end

function output = evaluateWell(cvId)
    global wel
    output = zeros(1,2);
    wellId = wel.cv2well(cvId);
    if ~isnan(wellId)
        output(1) = wel.qo(wellId);
        output(2) = wel.qw(wellId);
    end
end

function output = vectorizedWell()
    global wel msh
    output = zeros(msh.nMesh,2);
    output(wel.well2cv,:) = [wel.qo;wel.qw]';
end