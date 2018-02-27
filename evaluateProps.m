function [ outMesh ] = evaluateProps(mesh)
%EVALUATEPROPS Summary of this function goes here
%   Detailed explanation goes here
    
    outMesh = upwindConnections(mesh);

    satFunc = saturationFunctions;
    pvtFunc = pvtFunctions;
    compFunc = compactionFunctions;
    
    phase=phases.OIL;
    outMesh.bo  = pvtFunc.fvf(mesh.p,phase,mesh.PVTNUM);
    outMesh.muo = pvtFunc.viscosity(mesh.p,phase,mesh.PVTNUM);
    outMesh.roo = pvtFunc.density(mesh.p,phase,mesh.PVTNUM);
    outMesh.kro = satFunc.kr(mesh.so, mesh.sw,mesh.sg,phase,mesh.SATNUM);
    
    phase=phases.WATER;
    outMesh.bw  = pvtFunc.fvf(mesh.p,phase,mesh.PVTNUM);
    outMesh.muw = pvtFunc.viscosity(mesh.p,phase,mesh.PVTNUM);
    outMesh.row = pvtFunc.density(mesh.p,phase,mesh.PVTNUM);
    outMesh.krw = satFunc.kr(mesh.so, mesh.sw,mesh.sg,phase,mesh.SATNUM);
    outMesh.pcow = satFunc.pc(mesh.sw,phase,mesh.SATNUM);
    
    phase=phases.GAS;
    outMesh.bg  = pvtFunc.fvf(mesh.p,phase,mesh.PVTNUM);
    outMesh.mug = pvtFunc.viscosity(mesh.p,phase,mesh.PVTNUM);
    outMesh.rog = pvtFunc.density(mesh.p,phase,mesh.PVTNUM);
    outMesh.krg = satFunc.kr(mesh.so, mesh.sw,mesh.sg,phase,mesh.SATNUM);
    outMesh.pcog = satFunc.pc(mesh.sg,phase,mesh.SATNUM);
    
    outMesh.poreVolume = outMesh.poreVolume0 .* compFunc.porvMult(mesh.p,mesh.PVTNUM);
end

function [outMesh] = upwindConnections(mesh)
    outMesh = mesh;
    outMesh.upWind = mesh.face2cv(:,1).*(mesh.p(mesh.face2cv(:,1))>=mesh.p(mesh.face2cv(:,2)))'+ mesh.face2cv(:,2).*(mesh.p(mesh.face2cv(:,1))<mesh.p(mesh.face2cv(:,2)))';
end
