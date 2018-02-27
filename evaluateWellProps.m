function [ wellOut ] = evaluateWellProps( wells, mesh)
    wellOut = wells;
    for wellId=1:wells.nWells
        WI = wells.wellIndex(wellId);
        wellCv = wells.well2cv(wellId); po = mesh.p(wellCv);
        kro = mesh.kro(wellCv); krw = mesh.krw(wellCv);
        bo = mesh.bo(wellCv); bw = mesh.bw(wellCv);
        muo = mesh.muo(wellCv); muw = mesh.muw(wellCv);
        if wells.wellType(wellId) == wellTypes.PRODUCER
            oilMobility = kro/bo/muo;
            waterMobility = krw/bw/muw;
        else
            oilMobility = kro/bo/muo+krw/bw/muw;
            waterMobility = oilMobility;
        end
        if wells.activeConstraint(wellId)==wellConstraints.BHP
            pbh = wells.constraintValue(wellId);
            wellOut.pbh(wellId) = pbh;
            if wells.wellType(wellId) == wellTypes.PRODUCER
                wellOut.qo(wellId) = WI*oilMobility*(pbh-po);
                wellOut.qw(wellId) = WI*waterMobility*(pbh-po);
            else
                if wells.wellPhase(wellId) == phases.OIL
                    wellOut.qo(wellId) = WI*oilMobility*(pbh-po);
                    wellOut.qw(wellId) = 0.0;
                elseif wells.wellPhase(wellId) == phases.WATER
                    wellOut.qo(wellId) = 0.0;
                    wellOut.qw(wellId) = WI*waterMobility*(pbh-po);
                end
            end
        elseif wells.activeConstraint(wellId)==wellConstraints.ORAT
            qo = wells.constraintValue(wellId);
            wellOut.qo(wellId) = qo;
            pbh = qo/(WI*oilMobility) + po;
            wellOut.pbh(wellId) = pbh;
            if wells.wellType(wellId)==wellTypes.PRODUCER
                wellOut.qw(wellId) = WI*waterMobility*(pbh-po);
            else
                wellOut.qw(wellId) = 0.0;
            end
        elseif wells.activeConstraint(wellId)==wellConstraints.WRAT
            qw = wells.constraintValue(wellId);
            wellOut.qw(wellId) = qw;
            pbh = qw/(WI*waterMobility)+po;
            wellOut.pbh(wellId) = pbh;
            if wells.wellType(wellId)==wellTypes.PRODUCER
                wellOut.qo(wellId) = WI*oilMobility*(pbh-po);
            else
                wellOut.qo(wellId) = 0.0;
            end
        end
    end
    
end

