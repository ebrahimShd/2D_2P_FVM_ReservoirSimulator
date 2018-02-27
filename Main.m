clear;
tic
unit = units.FIELD;
dt = 0.0001*toSI(unit,measurements.time); 
endTime=2*toSI(unit,measurements.time); maxItr=20;


boundary = [0,0;1000,0;1000,10;0 10; 0,0] .* toSI(unit, measurements.distance);
nX = 1000; nY = 1;nZ = 1; nMesh = nX*nY*nZ;
dX = 1.*ones(1,nX).* toSI(unit, measurements.distance);
dY = 10.*ones(1,nY).* toSI(unit, measurements.distance);
dZ = 10.* toSI(unit, measurements.distance);

porosity = 0.2*ones(1,nMesh);
permX = 2000*ones(1,nMesh); permX = permX.*toSI(unit,measurements.permeability);
permY = permX; permZ = 0.1.*permX;
permeability = perm2tensor(permX,permY,permZ);
NTG = 1*ones(1,nMesh);

p0  = 4000*ones(1,nMesh).*toSI(unit,measurements.pressure);
sw0 = 0.2*ones(1,nMesh);
sg0 = zeros(1,nMesh);
nUpwind = 1;
upwindFormula = @singlePointUpwind;
nPhases = 2;
pMaxChangeTol = 100*toSI(unit,measurements.pressure);
sMaxChangeTol = 0.2;
easyItre = 9;
maxTS = 100*toSI(unit,measurements.time);

noWells = 2;
wellPosition=[50,5;950,5].*toSI(unit,measurements.distance);
wellRw=[0.26,0.26].*toSI(unit,measurements.distance);
wellType = [wellTypes.INJECTOR,wellTypes.PRODUCER];
wellPhase = [phases.WATER,phases.OIL];
wellConstraint = [wellConstraints.WRAT,wellConstraints.BHP];
wellConstraintValue = [5500*toSI(unit,measurements.liquidRate),2000*toSI(unit,measurements.pressure)];
wellIndex = [3,3].*toSI(unit,measurements.transmissibility);

geoMesh = meshGenerator2D(boundary,[nX,nY],dZ,dX,dY); %plot2Dmesh(geoMesh);
dynaMesh = meshPreProcess(geoMesh,porosity,permeability,NTG);
dynaMesh = initializeMesh(dynaMesh,false,p0,sw0,sg0,nUpwind,upwindFormula,nPhases,pMaxChangeTol,sMaxChangeTol);
wells = initializeWells(noWells,wellPosition,geoMesh,dynaMesh,wellRw,wellType,wellPhase,wellConstraint,wellConstraintValue,wellIndex);

time = 0;timeStep = 1;converged = false;
dynaMesh = evaluateProps(dynaMesh);
dynaMesh = gotoNextStep(dynaMesh,dt);
wells = evaluateWellProps(wells,dynaMesh);
dynaMesh = evaluateDerivatives(dynaMesh,wells,true);
residuals = evaluateResiduals(dynaMesh,wells,true);
results = initializeResults(dynaMesh, wells,unit);
outBoundArray = zeros(nPhases*nMesh,1);
chopped = 0;
while time<=endTime
    for itr = 1:maxItr
        if itr>1 || time==0 ,jacobianMatrix = assembleJacobi(dynaMesh);end
        %plotJacobi(jacobianMatrix,1);
        residuals = -1.*residuals;
        [jacobianMatrix,residuals] = applyConstraint(dynaMesh, jacobianMatrix , residuals,outBoundArray);
        rcondJacobi = 1/condest(jacobianMatrix);
        if abs(rcondJacobi) < eps
            fprintf('The rcond jacobi after apply constraint: %d. Linear non-convergency.\n',rcondJacobi);
            break;
        end
        deltaSol = jacobianMatrix\residuals;
        [dynaMesh, isOut, outBoundArray] = loadResult(dynaMesh,deltaSol);
        if isOut, fprintf('@--a solution is out of bound in itr: %d\n',itr);end
        %plotResults3D(geoMesh,dynaMesh,nX,nY,unit);
        %pause(0.1);
        dynaMesh = evaluateProps(dynaMesh);
        wells = evaluateWellProps(wells,dynaMesh);
        dynaMesh = evaluateDerivatives(dynaMesh,wells,true);
        residuals = evaluateResiduals(dynaMesh,wells,true);
        converged = checkConvergency(residuals);
        %fprintf('itr %d ',itr);
        if converged
            time = time + dt;
            results = evaluateResults(results, dynaMesh, wells,unit,time);
            fprintf('The time step %d is converged in %d iteration with dt: %f days. Time: %f days\n',timeStep, itr, dt/toSI(unit,measurements.time), time/toSI(unit,measurements.time));
%             plotResults3D(geoMesh,dynaMesh,nX,nY,unit);
            plotResults(results);
            pause(0.01);
            break;
        end
    end
%     fprintf('**** dt = %f   easyItre = %d    chopped = %f \n',dt/toSI(unit,measurements.time), easyItre, chopped);
    if converged
        if itr<easyItre 
            dt = min(dt*(2-0.9*chopped),maxTS);
            if chopped>0, easyItre = min(easyItre*1.1, 8);end
            chopped = max(0,chopped-0.1);
        end
        dynaMesh = gotoNextStep(dynaMesh,dt);
        timeStep = timeStep + 1;
    else
        dt = dt/2;
        chopped = 1;
        easyItre = max(easyItre/2,5);
        fprintf('@-- Time step %d is back tracked. Time step is chopped to %f days.\n',timeStep, dt/toSI(unit,measurements.time));
        dynaMesh = backTrack(dynaMesh,dt);
        dynaMesh = evaluateProps(dynaMesh);
        wells = evaluateWellProps(wells,dynaMesh);
        dynaMesh = evaluateDerivatives(dynaMesh,wells,true);
        residuals = evaluateResiduals(dynaMesh,wells,true);
    end
end
toc


