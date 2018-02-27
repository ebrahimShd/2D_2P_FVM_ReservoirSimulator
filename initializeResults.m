function results = initializeResults(mesh,wells,unit)
	nWells = wells.nWells;
	results = struct('time',[],'wellQs', [], 'wellPbh', [], 'cvPo', [], 'cvSw',[]);
	wellQs = cell(1,nWells);
	wellPbh = cell(1,nWells);
	cvPo = mesh.p ./toSI(unit,measurements.pressure);
	cvSw = mesh.sw;
    for wellId = 1:nWells
        wellQs{wellId} = [wells.qo(wellId),wells.qw(wellId)]./toSI(unit,measurements.liquidRate);
        wellPbh{wellId} = wells.pbh(wellId)/toSI(unit,measurements.pressure);
    end
    results.time = 0;
    results.wellQs = wellQs;
    results.wellPbh = wellPbh;
    results.cvPo = cvPo;
    results.cvSw = cvSw;
end