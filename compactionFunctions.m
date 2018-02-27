classdef compactionFunctions
    %COMPACTIONFUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cf = [2000*toSI(units.FIELD,measurements.pressure), 5.6e-6*toSI(units.FIELD,measurements.compressibility)];
    end
    
    methods(Access=public)
        
        function out = porvMult(obj,pressure,region)
            out = exp(obj.cf(region,2)'.*(pressure-obj.cf(region,1)'));
        end
        
        function out = dPorvMultdP(obj,pressure,region)
            out = obj.cf(region,2).*exp(obj.cf(region,2).*(pressure'-obj.cf(region,1)));
            out = out';
        end
        
    end
    
end

