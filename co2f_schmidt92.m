function [sch,sch_valid] = co2f_schmidt92(temp,gas)
%co2f_schmidt92 Calculate Schmidt number
% Part of co2flux() toolbox - https://github.com/mvdh7/co2flux
% Source: Wanninkhof, 1992, J Geophys Res 97. doi:10.1029/92JC00188
% Inputs: <temp> = temperature / degC; <gas> = which gas.
% <gas> options (from Table A1):
% 1 or 'co2' = carbon dioxide CO2
% 2 or 'o2'  = oxygen O2
% Written by Matthew P. Humphreys, last updated 2016-11-07

switch lower(gas)
    case {1 'co2'}  % CO2, seawater salinity = 35
        A = 2073.1; B = 125.62; C = 3.6276; D = 0.043219;
    case {2 'o2'} % O2, seawater salinity = 35
        A = 1953.4; B = 128.00; C = 3.9918; D = 0.050091;
end %switch

sch = A - B*temp + C*temp.^2 - D*temp.^3;
sch_valid = temp > 0 & temp < 30;
end %function co2f_schmidt92