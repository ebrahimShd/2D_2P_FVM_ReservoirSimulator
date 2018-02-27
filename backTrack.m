function outMesh = backTrack(mesh,dt)
    outMesh = mesh;
    outMesh.dt = dt;
    outMesh.p = mesh.pPrev;
    outMesh.sw = mesh.swPrev;
    outMesh.so = mesh.soPrev;
    outMesh.sg = 1.- outMesh.sw - outMesh.so;
end
