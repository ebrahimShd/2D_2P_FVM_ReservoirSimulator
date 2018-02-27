function [jMat, res] = applyConstraint(mesh, jMatrix, residual, outBoundArray)
    
    jMat = jMatrix;
	res = residual;
	nCv = mesh.nMesh;
    nPhases = mesh.nPhases;
	jMat = jMat .* sparse((1.-outBoundArray)*ones(1,nCv*nPhases));
	jMat = jMat + diag(outBoundArray);
    res = res.*(1.-outBoundArray);
end
