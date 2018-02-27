function [ outmesh ] = initializeMesh( mesh,equilibrate,varargin)
    global g;
    g = 9.8;
    outmesh=mesh;
    outmesh.p=[];
    outmesh.sw=[];
    outmesh.sg=[];
    outmesh.so=[];
    if~equilibrate
        p = cell2mat(varargin(1));
        sw = cell2mat(varargin(2));
        sg = cell2mat(varargin(3));
        outmesh.nUpwind = cell2mat(varargin(4));
        outmesh.upwindFormula = cell2mat(varargin(5));
        outmesh.nPhases = cell2mat(varargin(6));
        outmesh.pMaxChangeTol = cell2mat(varargin(7));
        outmesh.sMaxChangeTol = cell2mat(varargin(8));
        outmesh.p=p;
        outmesh.sw=sw;
        outmesh.sg=sg;
        outmesh.so=(1.-(sw+sg));
    else
        tableDimen = 100;g = 9.8;maxItr=10;
        datum = cell2mat(varargin(1));pref  = cell2mat(varargin(2));woc=cell2mat(varargin(3));goc=cell2mat(varargin(4));
        outmesh.nUpwind = cell2mat(varargin(5)); 
        outmesh.upwindFormula = cell2mat(varargin(6));
        outmesh.nPhases = cell2mat(varargin(7));
        outmesh.pMaxChangeTol = cell2mat(varargin(8));
        outmesh.sMaxChangeTol = cell2mat(varargin(9));
        minDepth = min(mesh.depth);maxDepth = max(mesh.depth);
        nRegion = max(mesh.PVTNUM);
        
        minDepth = min(minDepth, min([datum,woc,goc]));
        maxDepth = max(maxDepth, max([datum,woc,goc]));
        depthTable =(minDepth:(maxDepth-minDepth)/(tableDimen-1):maxDepth);
        [depthTable,id]=putValue2Vactor(depthTable,goc);
        [depthTable,id]=putValue2Vactor(depthTable,woc);
        [depthTable,id]=putValue2Vactor(depthTable,datum);
        idDatum = id; idWoc = find(depthTable==woc,1); idGoc = find(depthTable==goc,1);
        
        pressureTable = cell(1,nRegion);
        for i=1:nRegion
            pressureTable{i} = propagatePressureFirst(depthTable, idDatum, idWoc, idGoc, pref, i);
        end
        for itr=1:maxItr
            for i=1:nRegion
                pressureTable{i} = propagatePressure(depthTable, pressureTable{i}, idDatum, idWoc, idGoc, pref, i);
            end
        end
        
        sf = saturationFunctions;
        for i = 1:mesh.nMesh
            region = mesh.PVTNUM(i);
            depth = mesh.depth(i);
            SWL = sf.getEndPoint(region,'SWL');
            SWU = sf.getEndPoint(region,'SWU');
            SGL = sf.getEndPoint(region,'SGL');
            SGU = sf.getEndPoint(region,'SGU');
            po = interp1q(depthTable',pressureTable{region}(:,1),depth);
            pw = interp1q(depthTable',pressureTable{region}(:,2),depth);
            pg = interp1q(depthTable',pressureTable{region}(:,3),depth);
            if po>pw
                sw = SWL;
                if pg>po,sg = SGU;else sg = SGL;end
            else
                sw=SWU;sg=SGL;
            end
            so=1-sw-sg;
            outmesh.p(i)=po;
            outmesh.sw(i)=sw;
            outmesh.sg(i)=sg;
            outmesh.so(i)=so;
        end
    end
end

function [arr,id]=putValue2Vactor(array,value)
    id = 0;idUp=0;arrLength = length(array);
    for i=1:arrLength,if array(i)==value,id = i;break;elseif array(i)>value,idUp=i;break;end;end
    if id>0,arr=array;return;end
    if idUp<=1,id=[];arr=[];return;end
    
    id = idUp-1;arr = NaN(1,arrLength+1);
    for i=1:arrLength+1
        if i<idUp
            arr(i) = array(i);
        elseif i==idUp
            arr(i)=value;
        else
            arr(i)=array(i-1);
        end
    end
end

function press=propagatePressureFirst(depthTable, idDatum, idWoc, idGoc, pref, region)
    global g;
    lenTable = length(depthTable);
    po = NaN(lenTable,1);pw = NaN(lenTable,1);pg = NaN(lenTable,1);
    pvtFunc=pvtFunctions;
    %oil pressure datum down
    for i=idDatum:-1:1
        if i == idDatum
            po(i)=pref;
        else
            po(i) = po(i+1)+pvtFunc.density(po(i+1),'OIL',region)*g*(depthTable(i+1)-depthTable(i));
        end
    end
    %oil pressure datum up
    for i=idDatum+1:lenTable
        po(i) = po(i-1)-pvtFunc.density(po(i-1),'OIL',region)*g*(depthTable(i)-depthTable(i-1));
    end
    
    %water pressure woc down
    for i=idWoc:-1:1
        if i == idWoc
            pw(i)=po(i);
        else
            pw(i) = pw(i+1)+pvtFunc.density(pw(i+1),'WATER',region)*g*(depthTable(i+1)-depthTable(i));
        end
    end
    %water pressure woc up
    for i=idWoc+1:lenTable
        pw(i) = pw(i-1)-pvtFunc.density(pw(i-1),'WATER',region)*g*(depthTable(i)-depthTable(i-1));
    end
    %gas pressure goc up
    for i=idGoc:lenTable
        if i == idGoc
            pg(i)=po(i);
        else
            pg(i) = pg(i-1)-pvtFunc.density(pg(i-1),'GAS',region)*g*(depthTable(i)-depthTable(i-1));
        end
    end
    %gas pressure goc down
    for i=idGoc-1:-1:1
        pg(i) = pg(i+1)+pvtFunc.density(pg(i+1),'GAS',region)*g*(depthTable(i+1)-depthTable(i));
    end
    press = [po,pw,pg];
end

function press=propagatePressure(depthTable, pressureTable, idDatum, idWoc, idGoc, pref, region)
    global g;
    lenTable = length(depthTable);
    po = pressureTable(:,1);pw = pressureTable(:,2);pg = pressureTable(:,3);
    pvtFunc = pvtFunctions;
    
    %oil pressure datum down
    for i=idDatum:-1:1
        if i == idDatum
            po(i)=pref;
        else
            den=[pvtFunc.density(po(i),'OIL',region),pvtFunc.density(po(i+1),'OIL',region)];
            po(i) = po(i+1)+mean(den)*g*(depthTable(i+1)-depthTable(i));
        end
    end
    %oil pressure datum up
    for i=idDatum+1:lenTable
        den=[pvtFunc.density(po(i-1),'OIL',region),pvtFunc.density(po(i),'OIL',region)];
        po(i) = po(i-1)-mean(den)*g*(depthTable(i)-depthTable(i-1));
    end
    
    %water pressure woc down
    for i=idWoc:-1:1
        if i == idWoc
            pw(i)=po(i);
        else
            den=[pvtFunc.density(pw(i),'WATER',region),pvtFunc.density(pw(i+1),'WATER',region)];
            pw(i) = pw(i+1)+mean(den)*g*(depthTable(i+1)-depthTable(i));
        end
    end
    %water pressure woc up
    for i=idWoc+1:lenTable
        den=[pvtFunc.density(pw(i-1),'WATER',region),pvtFunc.density(pw(i),'WATER',region)];
        pw(i) = pw(i-1)-mean(den)*g*(depthTable(i)-depthTable(i-1));
    end
    %gas pressure goc up
    for i=idGoc:lenTable
        if i == idGoc
            pg(i)=po(i);
        else
            den = [pvtFunc.density(pg(i-1),'GAS',region),pvtFunc.density(pg(i),'GAS',region)];
            pg(i) = pg(i-1)-mean(den)*g*(depthTable(i)-depthTable(i-1));
        end
    end
    %gas pressure goc down
    for i=idGoc-1:-1:1
        den = [pvtFunc.density(pg(i),'GAS',region),pvtFunc.density(pg(i+1),'GAS',region)];
        pg(i) = pg(i+1)+mean(den)*g*(depthTable(i+1)-depthTable(i));
    end
    press = [po,pw,pg];
end

