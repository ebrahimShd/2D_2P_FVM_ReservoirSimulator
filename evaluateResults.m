function outResults = evaluateResults(results,mesh,wells,unit,time)
	outResults = results;
	nMesh = mesh.nMesh;
	nWells = wells.nWells;
	outResults.time = [results.time;time/toSI(unit,measurements.time)];
	for wellId = 1:nWells
		outResults.wellQs{wellId} = [results.wellQs{wellId};wells.qo(wellId)/toSI(unit,measurements.liquidRate),wells.qw(wellId)/toSI(unit,measurements.liquidRate)];
		outResults.wellPbh{wellId} = [results.wellPbh{wellId};wells.pbh(wellId)/toSI(unit,measurements.pressure)];;
	end
	outResults.cvPo = [results.cvPo;mesh.p ./toSI(unit,measurements.pressure)];
	outResults.cvSw = [results.cvSw;mesh.sw];
end