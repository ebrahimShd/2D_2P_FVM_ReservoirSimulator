function [ outMesh ] = evaluateDerivatives( mesh,wells,vectorized )
    if vectorized
        outMesh = vectorizedTempDerivatives(mesh);
        outMesh = vectorizedSpacDerivatives(outMesh);
        outMesh = vectorizedWellDerivatives(outMesh,wells);
    else
        outMesh = evaluateTempDerivatives(mesh);
        outMesh = evaluateSpacDerivatives(outMesh);
        outMesh = evaluateWellDerivatives(outMesh,wells);
    end
end

function outMesh = evaluateTempDerivatives(mesh)
    outMesh=mesh;
    der = cell(1,4); for i=1:4, der{i} = zeros(mesh.nMesh,1);end
    dt = mesh.dt;
    for cvId=1:mesh.nMesh
        compFunc = compactionFunctions;
        pvtFunc = pvtFunctions;
        satFunc = saturationFunctions;
        pRegion = mesh.PVTNUM(cvId);
        sRegion = mesh.SATNUM(cvId);
        porv0=mesh.poreVolume0(cvId);
        porv=mesh.poreVolume(cvId);
        po=mesh.p(cvId);
        pw=po-mesh.pcow(cvId);
        sw = mesh.sw(cvId);
        so = mesh.so(cvId);
        bw = mesh.bw(cvId);
        bo = mesh.bo(cvId);
        
        dPvdP=porv0*compFunc.dPorvMultdP(po,pRegion);
        dBodP=pvtFunc.dFvfdP(po,phases.OIL,pRegion);
        dBwdP=pvtFunc.dFvfdP(pw,phases.WATER,pRegion);
        dPcowdSw=satFunc.dPcdSw(sw,phases.WATER,sRegion);
        dBwdSw=dBwdP*(-dPcowdSw);
        
        dOtempdP = so*(dPvdP*bo-dBodP*porv)/(bo^2)/dt;
        dWtempdP = sw*(dPvdP*bw-dBwdP*porv)/(bw^2)/dt;
        dOtempdSw = -1*porv/bo/dt;
        dWtempdSw = porv*(bw - sw*dBwdSw)/(bw^2)/dt;
        
        der{1}(cvId)=dOtempdP;der{2}(cvId)=dOtempdSw;
        der{3}(cvId)=dWtempdP;der{4}(cvId)=dWtempdSw;
    end
    outMesh.temporalDerivatives = der;
end

function outMesh = evaluateSpacDerivatives(mesh)
%     [outMesh]=upwindConnections(mesh);
    outMesh = mesh;
    [dTrd,mobility] = getdTrd(outMesh);
    [dDpd,dp] = getdDpd(outMesh);
    outMesh.mobility = mobility;
    outMesh.deltaP = dp;
    outMesh.dTrd = dTrd;
    outMesh.dDpd = dDpd;
end

function [der,mobility]=getdTrd(mesh)
    nUpwind = mesh.nUpwind; upWind = mesh.upWind; nFace=mesh.nFace;
    der = cell(1,4); for i=1:4,der{i} = zeros(mesh.nFace,nUpwind);end
    mobility = zeros(nFace, 2);
    pvtFunc = pvtFunctions;                      
    satFunc = saturationFunctions;
    oilMobility = zeros(1,nUpwind);     waterMobility = zeros(1,nUpwind);
    for faceId=1:mesh.nFace
        faceCvs=upWind(faceId,:);
        for cvCounter=1:nUpwind
            poUp = mesh.p(faceCvs(cvCounter));           pwUp = poUp - mesh.pcow(faceCvs(cvCounter));
            soUp = mesh.so(faceCvs(cvCounter));          swUp = mesh.sw(faceCvs(cvCounter));            sgUp = mesh.sg(faceCvs(cvCounter));
            pRegionUp = mesh.PVTNUM(faceCvs(cvCounter)); sRegionUp = mesh.SATNUM(faceCvs(cvCounter));
            kroUp = mesh.kro(faceCvs(cvCounter));        krwUp = mesh.krw(faceCvs(cvCounter));
            boUp = mesh.bo(faceCvs(cvCounter));          bwUp = mesh.bw(faceCvs(cvCounter));
            muoUp = mesh.muo(faceCvs(cvCounter));        muwUp = mesh.muw(faceCvs(cvCounter));
            oilMobility(cvCounter) = kroUp/muoUp/boUp;   waterMobility(cvCounter) = krwUp/muwUp/bwUp;
            
            dMuodP = pvtFunc.dViscdP(poUp,phases.OIL, pRegionUp);            dMuwdP = pvtFunc.dViscdP(pwUp,phases.WATER, pRegionUp);
            dBodP = pvtFunc.dFvfdP(poUp,phases.OIL, pRegionUp);              dBwdP = pvtFunc.dFvfdP(pwUp,phases.WATER, pRegionUp);
            dKrodSw = satFunc.dKrdSw(soUp,swUp,sgUp, phases.OIL,sRegionUp);  dKrwdSw = satFunc.dKrdSw(soUp,swUp,sgUp, phases.WATER,sRegionUp);
            dMuwdSw = -dMuwdP*satFunc.dPcdSw(swUp, phases.WATER,sRegionUp);  dBwdSw = -dBwdP*satFunc.dPcdSw(swUp, phases.WATER,sRegionUp);
            
            dOmobilitydP = -kroUp*(dMuodP*boUp+dBodP*muoUp)/(muoUp*boUp)^2; dOmobilitydSw = dKrodSw*(1/muoUp/boUp);
            dWmobilitydP = -krwUp*(dMuwdP*bwUp+dBwdP*muwUp)/(muwUp*bwUp)^2; dWmobilitydSw = (dKrwdSw*muwUp*bwUp-krwUp*(dMuwdSw*bwUp+dBwdSw*muwUp))/(muwUp*bwUp)^2;
            
            der{1}(faceId,cvCounter)=dOmobilitydP;     der{2}(faceId,cvCounter)=dOmobilitydSw;
            der{3}(faceId,cvCounter)=dWmobilitydP;     der{4}(faceId,cvCounter)=dWmobilitydSw;
        end
        mO = mesh.upwindFormula(oilMobility);       mW = mesh.upwindFormula(waterMobility);
        mobility(faceId,1) = mO;                 mobility(faceId,2) = mW;
    end
end

function [der,dp]=getdDpd(mesh)
    der = cell(1,4);for i=1:4, der{i}=zeros(mesh.nFace,2);end
    pvtFunc = pvtFunctions;    satFunc = saturationFunctions;
    dp = zeros(mesh.nFace,2);
    for faceId = 1:mesh.nFace
        cv1 = mesh.face2cv(faceId,1);               cv2 = mesh.face2cv(faceId,2);
        po1 = mesh.p(cv1);                          po2 = mesh.p(cv2);
        pw1 = mesh.p(cv1)-mesh.pcow(cv1);           pw2 = mesh.p(cv2)-mesh.pcow(cv2);
        sw1 = mesh.sw(cv1);                         sw2 = mesh.sw(cv2);
        z1 = mesh.depth(cv1);                       z2 = mesh.depth(cv2);
        reg1 = mesh.PVTNUM(cv1);                    reg2 = mesh.PVTNUM(cv2);
        sreg1 = mesh.SATNUM(cv1);                   sreg2 = mesh.SATNUM(cv2);
        roO1 = pvtFunc.density(po1,phases.OIL,reg1);     roO2 = pvtFunc.density(po2,phases.OIL,reg2);
        roW1 = pvtFunc.density(pw1,phases.WATER,reg1);   roW2 = pvtFunc.density(pw2,phases.WATER,reg2);
        roOavr = roO1/2+roO2/2;                     roWavr = roW1/2+roW2/2;
        dRoodP1 = 0.5*pvtFunc.dDensdP(po1,phases.OIL,reg1);          dRoodP2 = 0.5*pvtFunc.dDensdP(po2,phases.OIL,reg2);
        dRowdP1 = 0.5*pvtFunc.dDensdP(pw1,phases.WATER,reg1);        dRowdP2 = 0.5*pvtFunc.dDensdP(pw2,phases.WATER,reg2);
        dPcowdSw1=satFunc.dPcdSw(sw1,phases.WATER,sreg1);            dPcowdSw2=satFunc.dPcdSw(sw2,phases.WATER,sreg2);
        
        dDpodP1 = -1+dRoodP1*(z2-z1);        dDpodSw1=  0.0;
        dDpodP2 =  1+dRoodP2*(z2-z1);        dDpodSw2=  0.0;
        dDpwdP1 = -1+dRowdP1*(z2-z1);        dDpwdSw1=  -dRowdP1*dPcowdSw1*(z2-z1)+dPcowdSw1;
        dDpwdP2 =  1+dRowdP2*(z2-z1);        dDpwdSw2=  -dRowdP2*dPcowdSw2*(z2-z1)-dPcowdSw2;
        
        der{1}(faceId,1)= dDpodP1;      der{1}(faceId,2)= dDpodP2;
        der{2}(faceId,1)= dDpodSw1;     der{2}(faceId,2)= dDpodSw2;
        der{3}(faceId,1)= dDpwdP1;      der{3}(faceId,2)= dDpwdP2;
        der{4}(faceId,1)= dDpwdSw1;     der{4}(faceId,2)= dDpwdSw2;
        
        dp(faceId,1) = (po2-po1)+roOavr*(z2-z1);
        dp(faceId,2) = (pw2-pw1)+roWavr*(z2-z1);
    end
end

function outMesh = evaluateWellDerivatives(mesh,wells)
    outMesh = mesh;
    der=cell(1,4); for i=1:4,der{i}=zeros(mesh.nMesh,1);end
	pvtFunc = pvtFunctions;         satFunc = saturationFunctions;
    for wellId=1:wells.nWells
        cvId = wells.well2cv(wellId);
        if wells.activeConstraint(wellId) == wellConstraints.BHP
            po = mesh.p(cvId);              pw = po - mesh.pcow(cvId);
            so = mesh.so(cvId);             sw = mesh.sw(cvId);         sg = mesh.sg(cvId);
            kro = mesh.kro(cvId);           krw = mesh.krw(cvId);
            bo = mesh.bo(cvId);             bw = mesh.bw(cvId);
            muo = mesh.muo(cvId);           muw = mesh.muw(cvId);
            pRegion = mesh.PVTNUM(cvId);    sRegion = mesh.SATNUM(cvId);
            Mo = kro/(bo*muo);              Mw = krw/(muw*bw);
            WI = wells.wellIndex(wellId);   pbh = wells.pbh(wellId);
            
            dMuodP = pvtFunc.dViscdP(po,phases.OIL, pRegion);
            dBodP = pvtFunc.dFvfdP(po,phases.OIL, pRegion);
            dMuwdP = pvtFunc.dViscdP(pw,phases.WATER, pRegion);
            dBwdP = pvtFunc.dFvfdP(pw,phases.WATER, pRegion);
            dKrodSw = satFunc.dKrdSw(so,sw,sg, phases.OIL,sRegion);
            dKrwdSw = satFunc.dKrdSw(so,sw,sg, phases.WATER,sRegion);
            dMuwdSw = -dMuwdP*satFunc.dPcdSw(sw, phases.WATER,sRegion);
            dBwdSw = -dBwdP*satFunc.dPcdSw(sw, phases.WATER,sRegion);
            
            if wells.wellType(wellId) == wellTypes.PRODUCER
                dModP = -kro*(dMuodP*bo+dBodP*muo)/(muo*bo)^2;
                dMwdP = -krw*(dMuwdP*bw+dBwdP*muw)/(muw*bw)^2;
                dModSw = dKrodSw*(1/muo/bo);
                dMwdSw = (dKrwdSw*muw*bw-krw*(dMuwdSw*bw+dBwdSw*muw))/(muw*bw)^2;
            else
                Mo = Mo+Mw; Mw = Mw+Mo;
                dModP  = -kro*(dMuodP*bo+dBodP*muo)/(muo*bo)^2 - krw*(dMuwdP*bw+dBwdP*muw)/(muw*bw)^2;
                dMwdP  = dModP;
                dModSw = dKrodSw*(1/muo/bo) + (dKrwdSw*muw*bw-krw*(dMuwdSw*bw+dBwdSw*muw))/(muw*bw)^2;
                dMwdSw = dModSw;
            end
            
            dWellodP  = WI*(dModP*(pbh-po)-Mo);     der{1}(cvId)=dWellodP;
            dWellodSw = WI*dModSw*(pbh-po);         der{2}(cvId)=dWellodSw;
            dWellwdP  = WI*(dMwdP*(pbh-po)-Mw);     der{3}(cvId)=dWellwdP;
            dWellwdSw = WI*dMwdSw*(pbh-po);         der{4}(cvId)=dWellwdSw;
            if wells.wellType(wellId) == wellTypes.INJECTOR && wells.wellPhase(wellId) == phases.OIL
                der{3}(cvId)=0.0;  der{4}(cvId)=0.0;
            end
            if wells.wellType(wellId) == wellTypes.INJECTOR && wells.wellPhase(wellId) == phases.WATER
                der{1}(cvId)=0.0;  der{2}(cvId)=0.0;
            end
        else
            dWellodP  = 0.0;     der{1}(cvId)=dWellodP;
            dWellodSw = 0.0;     der{2}(cvId)=dWellodSw;
            dWellwdP  = 0.0;     der{3}(cvId)=dWellwdP;
            dWellwdSw = 0.0;     der{4}(cvId)=dWellwdSw;
        end
        
        
    end
	outMesh.wellDerivatives = der;
end

function outMesh = vectorizedTempDerivatives(mesh)
    outMesh=mesh;
    %der = cell(1,4);
    dt = mesh.dt;
    compFunc = compactionFunctions;
    pvtFunc = pvtFunctions;
    satFunc = saturationFunctions;
    
    dPvdP=mesh.poreVolume0.*compFunc.dPorvMultdP(mesh.p,mesh.PVTNUM);
    dBodP=pvtFunc.dFvfdP(mesh.p,phases.OIL,mesh.PVTNUM);
    dBwdP=pvtFunc.dFvfdP(mesh.p,phases.WATER,mesh.PVTNUM);
    dPcowdSw=satFunc.dPcdSw(mesh.sw,phases.WATER,mesh.SATNUM);
    dBwdSw=dBwdP.*(-dPcowdSw);
        
    dOtempdP = mesh.so.*(dPvdP.*mesh.bo-dBodP.*mesh.poreVolume)./(mesh.bo.^2)./dt;
    dWtempdP = mesh.sw.*(dPvdP.*mesh.bw-dBwdP.*mesh.poreVolume)./(mesh.bw.^2)./dt;
    dOtempdSw = -1.*mesh.poreVolume./mesh.bo./dt;
    dWtempdSw = mesh.poreVolume.*(mesh.bw - mesh.sw.*dBwdSw)./(mesh.bw.^2)./dt;
    
%     der{1}=dOtempdP;der{2}=dOtempdSw;
%     der{3}=dWtempdP;der{4}=dWtempdSw;
    der = [dOtempdP;dOtempdSw;dWtempdP;dWtempdSw]';
    outMesh.temporalDerivatives = der;
end

function outMesh = vectorizedSpacDerivatives(mesh)
    outMesh = mesh;
    [dTrd,mobility] = vectorizeddTrd(outMesh);
    [dDpd,dp] = vectorizeddDpd(outMesh);
    outMesh.mobility = mobility;
    outMesh.deltaP = dp;
    outMesh.dTrd = dTrd;
    outMesh.dDpd = dDpd;
end

function [der,mobility]=vectorizeddTrd(mesh)
    nUpwind = mesh.nUpwind; 
    der = cell(1,4); for i=1:4, der{i}=zeros(mesh.nFace,nUpwind);end
    pvtFunc = pvtFunctions;                      
    satFunc = saturationFunctions;
    
    for i=1:nUpwind
        upWinds = mesh.upWind(:,i)';
        p = mesh.p(upWinds);
        sw = mesh.sw(upWinds);
        so = mesh.so(upWinds);
        sg = mesh.sg(upWinds);
        kro = mesh.kro(upWinds);
        krw = mesh.krw(upWinds);
        muo = mesh.muo(upWinds);
        muw = mesh.muw(upWinds);
        bo = mesh.bo(upWinds);
        bw = mesh.bw(upWinds);
        pw = p;
        pRegion = mesh.PVTNUM(mesh.upWind(i,:))';
        sRegion = mesh.SATNUM(mesh.upWind(i,:))';
        
        oilMobility = kro./muo./bo;                                   waterMobility = krw./muw./bw;
        dMuodP = pvtFunc.dViscdP(p,phases.OIL, pRegion);              dMuwdP = pvtFunc.dViscdP(pw,phases.WATER, pRegion);
        dBodP = pvtFunc.dFvfdP(p,phases.OIL, pRegion);                dBwdP = pvtFunc.dFvfdP(pw,phases.WATER, pRegion);
        dKrodSw = satFunc.dKrdSw(so,sw,sg, phases.OIL,sRegion);       dKrwdSw = satFunc.dKrdSw(so,sw,sg, phases.WATER,sRegion);
        dMuwdSw = -dMuwdP.*satFunc.dPcdSw(sw, phases.WATER,sRegion);  dBwdSw = -dBwdP.*satFunc.dPcdSw(sw, phases.WATER,sRegion);
        
        dOmobilitydP = -kro.*(dMuodP.*bo+dBodP.*muo)/(muo.*bo).^2; 
        dOmobilitydSw = dKrodSw.*(1./muo./bo);
        dWmobilitydP = -krw.*(dMuwdP.*bw+dBwdP.*muw)./(muw.*bw).^2;
        dWmobilitydSw = (dKrwdSw.*muw.*bw-krw.*(dMuwdSw.*bw+dBwdSw.*muw))./(muw.*bw).^2;
            
        der{1}(:,i)=dOmobilitydP';     der{2}(:,i)=dOmobilitydSw';
        der{3}(:,i)=dWmobilitydP';     der{4}(:,i)=dWmobilitydSw';
    end
        mO = mesh.upwindFormula(oilMobility);       mW = mesh.upwindFormula(waterMobility);
        mobility = [mO',mW'];
end

function [der,dp]=vectorizeddDpd(mesh)
    der = cell(1,4);for i=1:4, der{i}=zeros(mesh.nFace,2);end
    pvtFunc = pvtFunctions;    satFunc = saturationFunctions;
    dp = zeros(mesh.nFace,2);
    
    cv1 = mesh.face2cv(:,1)';                   cv2 = mesh.face2cv(:,2)';
    po1 = mesh.p(cv1);                          po2 = mesh.p(cv2);
    pw1 = mesh.p(cv1)-mesh.pcow(cv1);           pw2 = mesh.p(cv2)-mesh.pcow(cv2);
    sw1 = mesh.sw(cv1);                         sw2 = mesh.sw(cv2);
    z1 = mesh.depth(cv1)';                      z2 = mesh.depth(cv2)';
    reg1 = mesh.PVTNUM(cv1);                    reg2 = mesh.PVTNUM(cv2);
    sreg1 = mesh.SATNUM(cv1);                   sreg2 = mesh.SATNUM(cv2);
    roO1 = pvtFunc.density(po1,phases.OIL,reg1);     roO2 = pvtFunc.density(po2,phases.OIL,reg2);
    roW1 = pvtFunc.density(pw1,phases.WATER,reg1);   roW2 = pvtFunc.density(pw2,phases.WATER,reg2);
    roOavr = roO1./2+roO2./2;                     roWavr = roW1./2+roW2./2;
    dRoodP1 = 0.5*pvtFunc.dDensdP(po1,phases.OIL,reg1);          dRoodP2 = 0.5*pvtFunc.dDensdP(po2,phases.OIL,reg2);
    dRowdP1 = 0.5*pvtFunc.dDensdP(pw1,phases.WATER,reg1);        dRowdP2 = 0.5*pvtFunc.dDensdP(pw2,phases.WATER,reg2);
    dPcowdSw1=satFunc.dPcdSw(sw1,phases.WATER,sreg1);            dPcowdSw2=satFunc.dPcdSw(sw2,phases.WATER,sreg2);
    
    dDpodP1 = -1+dRoodP1.*(z2-z1);        dDpodSw1=  0.0;
    dDpodP2 =  1+dRoodP2.*(z2-z1);        dDpodSw2=  0.0;
    dDpwdP1 = -1+dRowdP1.*(z2-z1);        dDpwdSw1=  -dRowdP1.*dPcowdSw1.*(z2-z1)+dPcowdSw1;
    dDpwdP2 =  1+dRowdP2.*(z2-z1);        dDpwdSw2=  -dRowdP2.*dPcowdSw2.*(z2-z1)-dPcowdSw2;
    
    der{1}(:,1)= dDpodP1';      der{1}(:,2)= dDpodP2';
    der{2}(:,1)= dDpodSw1';     der{2}(:,2)= dDpodSw2';
    der{3}(:,1)= dDpwdP1';      der{3}(:,2)= dDpwdP2';
    der{4}(:,1)= dDpwdSw1';     der{4}(:,2)= dDpwdSw2';
    
    dp(:,1) = (po2-po1)+roOavr.*(z2-z1);
    dp(:,2) = (pw2-pw1)+roWavr.*(z2-z1);
end

function outMesh = vectorizedWellDerivatives(mesh,wells)
    outMesh = mesh;
    %der=cell(1,4); for i=1:4,der{i}=zeros(mesh.nMesh,1);end
    der = zeros(mesh.nMesh,4);
	pvtFunc = pvtFunctions;         satFunc = saturationFunctions;
    for wellId=1:wells.nWells
        cvId = wells.well2cv(wellId);
        if wells.activeConstraint(wellId) == wellConstraints.BHP
            po = mesh.p(cvId);              pw = po - mesh.pcow(cvId);
            so = mesh.so(cvId);             sw = mesh.sw(cvId);         sg = mesh.sg(cvId);
            kro = mesh.kro(cvId);           krw = mesh.krw(cvId);
            bo = mesh.bo(cvId);             bw = mesh.bw(cvId);
            muo = mesh.muo(cvId);           muw = mesh.muw(cvId);
            pRegion = mesh.PVTNUM(cvId);    sRegion = mesh.SATNUM(cvId);
            Mo = kro/(bo*muo);              Mw = krw/(muw*bw);
            WI = wells.wellIndex(wellId);   pbh = wells.pbh(wellId);
            
            dMuodP = pvtFunc.dViscdP(po,phases.OIL, pRegion);
            dBodP = pvtFunc.dFvfdP(po,phases.OIL, pRegion);
            dMuwdP = pvtFunc.dViscdP(pw,phases.WATER, pRegion);
            dBwdP = pvtFunc.dFvfdP(pw,phases.WATER, pRegion);
            dKrodSw = satFunc.dKrdSw(so,sw,sg, phases.OIL,sRegion);
            dKrwdSw = satFunc.dKrdSw(so,sw,sg, phases.WATER,sRegion);
            dMuwdSw = -dMuwdP*satFunc.dPcdSw(sw, phases.WATER,sRegion);
            dBwdSw = -dBwdP*satFunc.dPcdSw(sw, phases.WATER,sRegion);
            
            if wells.wellType(wellId) == wellTypes.PRODUCER
                dModP = -kro*(dMuodP*bo+dBodP*muo)/(muo*bo)^2;
                dMwdP = -krw*(dMuwdP*bw+dBwdP*muw)/(muw*bw)^2;
                dModSw = dKrodSw*(1/muo/bo);
                dMwdSw = (dKrwdSw*muw*bw-krw*(dMuwdSw*bw+dBwdSw*muw))/(muw*bw)^2;
            else
                Mo = Mo+Mw; Mw = Mw+Mo;
                dModP  = -kro*(dMuodP*bo+dBodP*muo)/(muo*bo)^2 - krw*(dMuwdP*bw+dBwdP*muw)/(muw*bw)^2;
                dMwdP  = dModP;
                dModSw = dKrodSw*(1/muo/bo) + (dKrwdSw*muw*bw-krw*(dMuwdSw*bw+dBwdSw*muw))/(muw*bw)^2;
                dMwdSw = dModSw;
            end
            
            dWellodP  = WI*(dModP*(pbh-po)-Mo);     der(cvId,1)=dWellodP; %der{1}(cvId)=dWellodP;
            dWellodSw = WI*dModSw*(pbh-po);         der(cvId,2)=dWellodSw;%der{2}(cvId)=dWellodSw;
            dWellwdP  = WI*(dMwdP*(pbh-po)-Mw);     der(cvId,3)=dWellwdP; %der{3}(cvId)=dWellwdP;
            dWellwdSw = WI*dMwdSw*(pbh-po);         der(cvId,4)=dWellwdSw;%der{4}(cvId)=dWellwdSw;
            if wells.wellType(wellId) == wellTypes.INJECTOR && wells.wellPhase(wellId) == phases.OIL
                %der{3}(cvId)=0.0;  der{4}(cvId)=0.0;
                der(cvId,3)=0.0;  der(cvId,4)=0.0;
            end
            if wells.wellType(wellId) == wellTypes.INJECTOR && wells.wellPhase(wellId) == phases.WATER
                %der{1}(cvId)=0.0;  der{2}(cvId)=0.0;
                der(cvId,1)=0.0;  der(cvId,2)=0.0;
            end
        else
            dWellodP  = 0.0;     der(cvId,1) = dWellodP; %der{1}(cvId)=dWellodP;
            dWellodSw = 0.0;     der(cvId,2) = dWellodSw;%der{2}(cvId)=dWellodSw;
            dWellwdP  = 0.0;     der(cvId,3) = dWellwdP; %der{3}(cvId)=dWellwdP;
            dWellwdSw = 0.0;     der(cvId,4) = dWellwdSw;%der{4}(cvId)=dWellwdSw;
        end
        
        
    end
	outMesh.wellDerivatives = der;
end









