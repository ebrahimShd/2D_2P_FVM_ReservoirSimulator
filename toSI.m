function value = toSI(from, measurement)
if from == units.SI
    value = 1;
    return
end

switch measurement
    case measurements.permeability
        switch from
            case units.METRIC, value = 9.869233e-16;
            case units.FIELD , value = 9.869233e-16;
        end
    case measurements.pressure
        switch from
            case units.METRIC, value = 100000;
            case units.FIELD , value = 6894.757;
        end
    case measurements.distance
        switch from
            case units.METRIC, value = 1;
            case units.FIELD , value = 0.3048;
        end
    case measurements.volume
        switch from
            case units.METRIC, value = 1;
            case units.FIELD , value = 0.1589873;
        end
    case measurements.gasFVF
        switch from
            case units.METRIC, value = 1;
            case units.FIELD , value = 5.614573/1000;
        end
    case measurements.viscosity
        switch from
            case units.METRIC, value = 0.001;
            case units.FIELD , value = 0.001;
        end
    case measurements.time
        switch from
            case units.METRIC, value = 86400;
            case units.FIELD , value = 86400;
        end
	case measurements.transmissibility
        switch from
            case units.METRIC, value = 0.001*1/100000/86400;
            case units.FIELD , value = 0.001*0.1589873/6894.757/86400;
        end
	case measurements.liquidRate
        switch from
            case units.METRIC, value = 1/86400;
            case units.FIELD , value = 0.1589873/86400;
        end
	case measurements.compressibility
        switch from
            case units.METRIC, value = 1/100000;
            case units.FIELD , value = 1/6894.757;
        end
	case measurements.density
        switch from
            case units.METRIC, value = 1;
            case units.FIELD , value = 0.45359237/0.028316846592;
        end
end
end


