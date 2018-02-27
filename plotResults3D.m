function []=plotResults3D(geoMesh,mesh,nX,nY,unit)
	p = mesh.p./toSI(unit,measurements.pressure);
	sw = mesh.sw; 
	cc = geoMesh.elementCenter;
	figure(2);
    colormap('jet');
	hold off;
	subplot(2,1,1);
	surf(reshape(cc(:,1),[nX,nY]),reshape(cc(:,2),[nX,nY]),reshape(p,[nX,nY]));
    view([-0.5 90]);
    grid('on');
    colorbar;
	title('pressure');
	subplot(2,1,2);
	surf(reshape(cc(:,1),[nX,nY]),reshape(cc(:,2),[nX,nY]),reshape(sw,[nX,nY]));
    view([-0.5 90]);
    grid('on');
    colorbar;
	title('water saturation');
    

end