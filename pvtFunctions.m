classdef pvtFunctions
    %PVTFUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %       Pref           FVF   C                   Vic  Viscosibility
        pvtw  = [2000*toSI(units.FIELD,measurements.pressure), 1.02, 							         1e-6*toSI(units.FIELD,measurements.compressibility), 0.05*toSI(units.FIELD,measurements.viscosity), 0*toSI(units.FIELD,measurements.compressibility) ];
        pvtg  = [2000*toSI(units.FIELD,measurements.pressure), 0.1*toSI(units.FIELD,measurements.gasFVF),1e-4*toSI(units.FIELD,measurements.compressibility), 0.01*toSI(units.FIELD,measurements.viscosity), 0*toSI(units.FIELD,measurements.compressibility) ];
        pvto  = [2000*toSI(units.FIELD,measurements.pressure), 1.3,  							         1e-5*toSI(units.FIELD,measurements.compressibility), 1.00*toSI(units.FIELD,measurements.viscosity), 0*toSI(units.FIELD,measurements.compressibility) ];
        rowsc = [71] *toSI(units.FIELD,measurements.density);
		rogsc = [0.1]*toSI(units.FIELD,measurements.density);
		roosc = [53] *toSI(units.FIELD,measurements.density);
    end
    
    methods(Access = public)
        
        function fvf = fvf( obj,pressure, phase, region )
            switch phase
                case phases.WATER
                    fvf = obj.getBw(pressure,region);
                case phases.GAS
                    fvf = obj.getBg(pressure,region);
                case phases.OIL
                    fvf = obj.getBo(pressure,region);
                otherwise
                    fvf = zeros(size(pressure));
                      
            end
        end
                
        function dens = density( obj,pressure, phase, region )
            fvfP = obj.fvf(pressure,phase,region);
            switch phase
                case phases.WATER
                    dens = obj.getRow(fvfP,region);
                case phases.GAS
                    dens = obj.getRog(fvfP,region);
                case phases.OIL
                    dens = obj.getRoo(fvfP,region);
                otherwise
                    dens = zeros(size(pressure));
            end
        end
        
        function visc = viscosity( obj, pressure, phase, region)
            switch phase
                case phases.WATER
                    visc = obj.getMuw(pressure,region);
                case phases.GAS
                    visc = obj.getMug(pressure,region);
                case phases.OIL
                    visc = obj.getMuo(pressure,region);
                otherwise
                    visc = zeros(size(pressure));
            end
        end
                
        function output = dFvfdP(obj, pressure, phase, region)
            switch phase
                case phases.WATER
                    output = obj.getdBwdP(pressure,region);
                case phases.GAS
                    output = obj.getdBgdP(pressure,region);
                case phases.OIL
                    output = obj.getdBodP(pressure,region);
                otherwise
                    output = zeros(size(pressure));
            end
        end
        
        function output = dViscdP(obj, pressure, phase, region)
            switch phase
                case phases.WATER
                    output = obj.getdMuwdP(pressure,region);
                case phases.GAS
                    output = obj.getdMugdP(pressure,region);
                case phases.OIL
                    output = obj.getdMuodP(pressure,region);
                otherwise
                    output = zeros(size(pressure));
            end
        end
        
        function output = dDensdP(obj, pressure, phase, region)
            switch phase
                case phases.WATER
                    output = obj.getdRowdP(pressure,region);
                case phases.GAS
                    output = obj.getdRogdP(pressure,region);
                case phases.OIL
                    output = obj.getdRoodP(pressure,region);
                otherwise
                    output = zeros(size(pressure));
            end
        end
    end
    
    methods(Access=private)
        
        function bw = getBw(obj,pw,region)
            bw = obj.pvtw(region,2)'.*exp(obj.pvtw(region,3)'.*(obj.pvtw(region,1)'-pw));
        end
        
        function bg = getBg(obj,pg,region)
            bg = obj.pvtg(region,2)'.*exp(obj.pvtg(region,3)'.*(obj.pvtg(region,1)'-pg));
        end
        
        function bo = getBo(obj,po,region)
            bo = obj.pvto(region,2)'.*exp(obj.pvto(region,3)'.*(obj.pvto(region,1)'-po));
        end
        
        function Row = getRow(obj,bw,region)
            Row = obj.rowsc(region,1)'./bw;
        end
        
        function Rog = getRog(obj,bg,region)
            Rog = obj.rogsc(region,1)'./bg;
        end
        
        function Roo = getRoo(obj,bo,region)
            Roo = obj.roosc(region,1)'./bo;
        end
        
        function muw = getMuw(obj,pw,region)
            muw = obj.pvtw(region,4)'.*exp(obj.pvtw(region,5)'.*(obj.pvtw(region,1)'-pw));
        end
        
        function mug = getMug(obj,pg,region)
            mug = obj.pvtg(region,4)'.*exp(obj.pvtg(region,5)'.*(obj.pvtg(region,1)'-pg));
        end
        
        function muo = getMuo(obj,po,region)
            muo = obj.pvto(region,4)'.*exp(obj.pvto(region,5)'.*(obj.pvto(region,1)'-po));
        end
        
        function out = getdBwdP(obj,pw,region)
            out = -obj.getBw(pw,region).*obj.pvtw(region,3)';
        end
        
        function out = getdBodP(obj,po,region)
            out = -obj.getBo(po,region).*obj.pvto(region,3)';
        end
        
        function out = getdBgdP(obj,pg,region)
            out = -obj.getBg(pg,region).*obj.pvtg(region,3)';
        end
        
        function out = getdMuwdP(obj,pw,region)
            out = -obj.getMuw(pw,region).*obj.pvtw(region,5)';
        end
        
        function out = getdMuodP(obj,po,region)
            out = -obj.getMuo(po,region).*obj.pvto(region,5)';
        end
        
        function out = getdMugdP(obj,pg,region)
            out = -obj.getMug(pg,region).*obj.pvtg(region,5)';
        end
        
        function out = getdRowdP(obj,pw,region)
            bw = obj.getBw(pw,region);
            dBwdP = obj.getdBwdP(pw,region);
            out = -obj.rowsc(region).*dBwdP./bw.^2;
        end
        
        function out = getdRoodP(obj,po,region)
            bo = obj.getBo(po,region);
            dBodP = obj.getdBodP(po,region);
            out = -obj.roosc(region).*dBodP./bo.^2;
        end
        
        function out = getdRogdP(obj,pg,region)
            bg = obj.getBg(pg,region);
            dBgdP = obj.getdBgdP(pg,region);
            out = -obj.rogsc(region).*dBgdP./bg.^2;
        end
        
        
    end
    
    
end

