classdef saturationFunctions
    
    properties (SetAccess = private)
        endPointTable = [0.2 0.201 0.17 1.0 0.0 0.01 0.0 0.8];
        coreyExpTable = [2 3 3 3];
    end
    
    methods (Access = public)
        
        function output = kr( obj, so, sw, sg, phase, region )
            switch phase
                case phases.WATER
                    output = obj.drainageKrw(sw,region);
                    output = output';
                case phases.GAS
                    output = obj.drainageKrg(sg,region);
                    output = output';
                case phases.OIL
                    krow = obj.drainageKrow(so,region);
                    krog = obj.drainageKrog(so,region);
                    SWL = obj.getEndPoint(region,endPoints.SWL);
                    SGL = obj.getEndPoint(region,endPoints.SGL);
                    sw = sw';sg = sg';
                    output1 = krow.*(sw<=SWL).*(sg<=SGL) + krog.*(sw<=SWL).*(sg>SGL) + krow.*(sw> SWL).*(sg<=SGL);
                    %output2 = ((sg.*krog + (sw-SWL).*krow)./(sg+sw-SWL)).*(sw>SWL).*(sg>SGL);
                    output = output1;% + output2;
                    output = output';
                otherwise
                    output = zeros(size(so));
            end
            
        end
        
        function output = pc( obj, saturation, phase, region )
            output = zeros(size(saturation));
        end
        
        function endPoint=getEndPoint(obj,region,EPN)
            endPoint = obj.endPointTable(region,EPN);
        end
        
        function out = dKrdSw(obj, so, sw, sg, phase, region)
            switch phase
                case phases.WATER
                    out = obj.dDrnKrwdSw(sw,region)';
                case phases.GAS
                    out = obj.dDrnKrgdSw(sg,region)';
                case phases.OIL
                    %krow = obj.drainageKrow(so,region)';
                    %krog = obj.drainageKrog(so,region)';
                    dKrowdSo = obj.dDrnKrowdSo(so,region);
                    SWL = obj.getEndPoint(region,endPoints.SWL);
                    SGL = obj.getEndPoint(region,endPoints.SGL);
                    SOWCR = obj.getEndPoint(region,endPoints.SOWCR);
                    sw = sw'; sg=sg';
                    out  = -1.*dKrowdSo.*(sw>SWL).*(sw<1-SOWCR-SGL).*(sg<=SGL);
                    %out2 = (((-1.*sw.*dKrowdSo+krow)-((sg-SGL).*krog + (sw-SWL).*krow))./(sg-SGL+sw-SWL).^2).*(sw>SWCR).*(sw<1-SOWCR-SGL).*(sg>SGL);
                    out = out';%+out2;
                otherwise
                    out = zeros(size(so));
            end
        end
        
        function out = dPcdSw(obj, saturation, phase, region)
            switch phase
                case phases.WATER
                    out = zeros(size(saturation));
                case phases.GAS
                    out = zeros(size(saturation));
                otherwise
                    out = zeros(size(saturation));
            end
        end
        
    end
    
    methods (Access = private)
        function [kr]=drainageKrw(obj,sw,region)
            SWCR = obj.getEndPoint(region,endPoints.SWCR);
            SOWCR = obj.getEndPoint(region,endPoints.SOWCR);
            SGL = obj.getEndPoint(region,endPoints.SGL);
            coreyExp = obj.getCoreyExponent(region,coreyPhases.WATER);
			maxKRW = 1.0;
            sw=sw';
            kr = maxKRW.*((sw-SWCR)./(1-SOWCR-SWCR-SGL)).^coreyExp.*(sw<=1.-SOWCR-SGL).*(sw>SWCR);
            kr = kr+maxKRW.*(sw>1.-SOWCR-SGL);
        end
        
        function [kr]=drainageKrg(obj, sg, region)
            SGCR = obj.getEndPoint(region,endPoints.SGCR);
            SOGCR = obj.getEndPoint(region,endPoints.SOGCR);
            SWL = obj.getEndPoint(region,endPoints.SWL);
            coreyExp = obj.getCoreyExponent(region,coreyPhases.GAS);
            maxKRG = 1.0;
            sg=sg';
            kr = maxKRG*((sg-SGCR)./(1-SOGCR-SGCR-SWL)).^coreyExp.*(sg<=1.-SOGCR-SWL).*(sg>SGCR);
            kr = kr+maxKRG.*(sg>1.-SOGCR-SWL);
        end
        
        function [kr]=drainageKrow(obj, so, region)
            SOWCR = obj.getEndPoint(region,endPoints.SOWCR);
            SGL = obj.getEndPoint(region,endPoints.SGL);
            SWL = obj.getEndPoint(region,endPoints.SWL);
            coreyExp = obj.getCoreyExponent(region,coreyPhases.OIL_WATER);
            maxKRO = 1.0;
            so=so';
            kr = maxKRO.*((so-SOWCR)./(1-SOWCR-SWL-SGL)).^coreyExp.*(so<=1.-SWL-SGL).*(so>SOWCR);
            kr = kr+maxKRO.*(so>1.-SWL-SGL);
        end
        
        function [kr]=drainageKrog(obj,so, region)
            SOGCR = obj.getEndPoint(region,endPoints.SOGCR);
            SGL = obj.getEndPoint(region,endPoints.SGL);
            SWL = obj.getEndPoint(region,endPoints.SWL);
            coreyExp = obj.getCoreyExponent(region,coreyPhases.OIL_GAS);
            maxKRO = 1.0;
            so=so';
            kr = maxKRO.*((so-SOGCR)./(1-SOGCR-SWL-SGL)).^coreyExp.*(so<=1.-SWL-SGL).*(so>SOGCR);
            kr = kr+maxKRO.*(so>1.-SWL-SGL);
        end
        
        function [coreyExp]=getCoreyExponent(obj,region,phase)
            coreyExp = obj.coreyExpTable(region,phase);
        end
        
        function out = dDrnKrwdSw(obj,sw,region)
            SWCR = obj.getEndPoint(region,endPoints.SWCR);
            SOWCR = obj.getEndPoint(region,endPoints.SOWCR);
            SGL = obj.getEndPoint(region,endPoints.SGL);
            coreyExp = obj.getCoreyExponent(region,coreyPhases.WATER);
			maxKRW = 1.0;sw = sw';
            out = maxKRW.*coreyExp.*(((sw-SWCR)./(1-SOWCR-SWCR-SGL)).^(coreyExp-1)).*(sw<=1-SOWCR-SGL).*(sw>SWCR);
        end
        
        function out = dDrnKrgdSw(obj,sw,region)
            out = zeros(size(sw))';
        end
        
        function out = dDrnKrowdSo(obj,so,region)
            SOWCR = obj.getEndPoint(region,endPoints.SOWCR);
            SGL = obj.getEndPoint(region,endPoints.SGL);
            SWL = obj.getEndPoint(region,endPoints.SWL);
            coreyExp = obj.getCoreyExponent(region,coreyPhases.OIL_WATER);
            maxKRO = 1.0;
            so = so';
            out = maxKRO.*coreyExp.*(((so-SOWCR)./(1-SOWCR-SWL-SGL)).^(coreyExp-1)).*(so<=1-SWL-SGL).*(so>SOWCR);
        end
    end
end

