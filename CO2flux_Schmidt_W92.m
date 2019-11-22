% CO2flux: air-sea CO2 fluxes with uncertainties.
% Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
function [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas)
%CO2flux_Schmidt_W92 Calculate Schmidt number following Wanninkhof (1992).
% Part of the CO2flux toolbox [https://github.com/mvdh7/CO2flux].
% Written by Matthew P. Humphreys [v1.0.0, last updated 2019-11-22].
% Inputs:
%   temp = Seawater temperature in degrees-C.
%   u10 = Wind speed at 10 m above the sea surface in m/s.
%   gas = Name of gas for Schmidt number calculation (case insensitive):
%     'CO2' = Carbon dioxide;
%     'O2' = Oxygen.
% Outputs:
%   Sch = The Schmidt number.
%   Sch_valid = Is the input temperature within the valid range?
switch lower(gas)
    case {1 'co2'}  % CO2, seawater salinity = 35
        A = 2073.1;
        B = 125.62;
        C = 3.6276;
        D = 0.043219;
    case {2 'o2'} % O2, seawater salinity = 35
        A = 1953.4;
        B = 128.00;
        C = 3.9918;
        D = 0.050091;
end
Sch = A - B*temp + C*temp.^2 - D*temp.^3;
Sch_valid = temp > 0 & temp < 30;
