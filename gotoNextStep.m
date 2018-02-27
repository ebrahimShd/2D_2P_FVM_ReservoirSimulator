function outMesh = gotoNextStep(mesh,dt)
    outMesh = mesh;
    outMesh.dt = dt;
    outMesh.pPrev = mesh.p;
    outMesh.swPrev = mesh.sw;
    outMesh.soPrev = mesh.so;
    outMesh.boPrev = mesh.bo;
    outMesh.bwPrev = mesh.bw;
    outMesh.pvPrev = mesh.poreVolume;
end