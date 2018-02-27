function [ jacobi ] = assembleMatrixNumerical( mesh,wells,originalResidual )
	varInc = [0.1,0.001];
	jacobi = zeros(2*mesh.nMesh,2*mesh.nMesh);
	for cvId=1:mesh.nMesh
		for var=1:2
			column = 2*(cvId-1)+var;
			if var ==1
				origP = mesh.p(cvId);
				mesh.p(cvId) = mesh.p(cvId)+varInc(var);
			else
				origSw = mesh.sw(cvId);
				mesh.sw(cvId) = mesh.sw(cvId)+varInc(var);
			end
			incResidual = calcResiduals(mesh,wells);
			if var == 1
				mesh.p(cvId)=origP;
			else
				mesh.sw(cvId)=origSw;
			end
			jacobi(:,column) = (incResidual-originalResidual)./varInc(var);	
		end
	end
	jacobi = sparse(jacobi);
end

function residuals=calcResiduals(mesh,well)
    global msh wel
    msh = mesh;wel = well;
    residuals = zeros(mesh.nMesh*2,1);
    for cvId = 1:mesh.nMesh
        temp = evaluateTemporal(cvId);
        spac = evaluateSpacial(cvId);
        well = evaluateWell(cvId);
        for phase = 1:2
            residuals(2*(cvId-1)+phase) = temp(phase)-spac(phase)-well(phase);
        end
    end
end

function output = evaluateTemporal(cvId)
    global msh
    output = zeros(1,2);
	pvtFunc = pvtFunctions;
	compFunc = compactionFunctions;
    sw = msh.sw(cvId);      				swPrev = msh.swPrev(cvId);
    so = 1-sw;								soPrev = msh.soPrev(cvId);
	po = msh.p(cvId);						pReg = msh.PVTNUM(cvId);
    bo = pvtFunc.fvf(po,'OIL',pReg);      	boPrev = msh.boPrev(cvId);
    bw = pvtFunc.fvf(po,'WATER',pReg);		bwPrev = msh.bwPrev(cvId);
	porvMult = compFunc.porvMult(po,pReg);
    pv = msh.poreVolume0(cvId)*porvMult;	pvPrev = msh.pvPrev(cvId);
    dt = msh.dt;
    output(1) = (pv*so/bo-pvPrev*soPrev/boPrev)/dt;
    output(2) = (pv*sw/bw-pvPrev*swPrev/bwPrev)/dt;
end

function output = evaluateSpacial(cvId)
    global msh
    output = zeros(1,2);
    cvFaces = msh.cv2face(cvId,:);cvFaces = cvFaces(~isnan(cvFaces)); nFace = length(cvFaces);
    for counter = 1:nFace
        faceId = cvFaces(counter);
        tran = msh.tran(faceId);
        faceCvs = msh.face2cv(faceId,:);
        mobility = getMobility(msh,faceId);
        dp = getDeltaP(msh,faceId); if cvId == faceCvs(2), dp = -1.*dp;end
        for phase=1:2
            output(phase) = output(phase)+tran*mobility(phase)*dp(phase);
        end
    end
end

function mobility = getMobility(mesh,faceId)
	satFunc = saturationFunctions; pvtFunc = pvtFunctions; 
	upwindCv = mesh.upWind(faceId);
    po = mesh.p(upwindCv);
	sw = mesh.sw(upwindCv);
	pReg = mesh.PVTNUM(upwindCv);
	sReg = mesh.SATNUM(upwindCv);
	so = 1-sw;
	sg = 0;
	kro = satFunc.kr(so,sw,sg,'OIL',sReg);
	krw = satFunc.kr(so,sw,sg,'WATER',sReg);
	bo = pvtFunc.fvf(po,'OIL',pReg);
	bw = pvtFunc.fvf(po,'WATER',pReg);
	muo = pvtFunc.viscosity(po,'OIL',pReg);
	muw = pvtFunc.viscosity(po,'WATER',pReg);
	mobility = [kro/bo/muo,krw/bw/muw];
end

function deltaP = getDeltaP(mesh,faceId)
	satFunc = saturationFunctions; pvtFunc = pvtFunctions; 
	faceCvs = mesh.face2cv(faceId,:);
	
	cv1 = faceCvs(1);		p1 = mesh.p(cv1);
	cv2 = faceCvs(2);		p2 = mesh.p(cv2);

	deltaP = [p2-p1,p2-p1];
end

function output = evaluateWell(cvId)
    global wel msh
    pvtFunc = pvtFunctions;
	satFunc = saturationFunctions;
	output = zeros(1,2);
    wellId = wel.cv2well(cvId);
	po = msh.p(cvId);
	sw = msh.sw(cvId);
	so = 1-sw;
	sg=0;
	sReg = msh.SATNUM(cvId);
	pReg = msh.PVTNUM(cvId);
	kro = satFunc.kr(so,sw,sg,'OIL',sReg);
	krw = satFunc.kr(so,sw,sg,'WATER',sReg);
	bo = pvtFunc.fvf(po,'OIL',pReg);
	bw = pvtFunc.fvf(po,'WATER',pReg);
	muo = pvtFunc.viscosity(po,'OIL',pReg);
	muw = pvtFunc.viscosity(po,'WATER',pReg);
    if ~isnan(wellId)
        WI = wel.wellIndex(wellId);
        if wel.activeConstraint(wellId)==wellConstraints.BHP
            pbh = wel.constraintValue(wellId);
            mobilityOil = kro/muo/bo;
            mobilityWater = krw/muw/bw;
            if wel.wellType(wellId) == wellTypes.PRODUCER
                output(1) = WI*mobilityOil*(pbh-po);
                output(2) = WI*mobilityWater*(pbh-po);
            else
                mobility = mobilityOil+mobilityWater;
                if wel.wellPhase(wellId)==phases.WATER
                    output(2) = WI*mobility*(pbh-po);
                    output(1) = 0.0;
                else
                    output(1) = WI*mobility*(pbh-po);
                    output(2) = 0.0;
                end
            end
        else
            output(1) = wel.qo(wellId);
            output(2) = wel.qw(wellId);
        end
    end
end







