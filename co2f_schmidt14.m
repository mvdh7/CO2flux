function [sch,sch_valid] = co2f_schmidt14(temp,gas)
%co2f_schmidt14 Calculate Schmidt number
% Part of co2flux() toolbox - https://github.com/mvdh7/co2flux
% Source: Wanninkhof, 2014, Limnol Oceanogr Methods 12, 351-362
% Inputs: <temp> = temperature / degC; <gas> = which gas.
% <gas> options (from 'seawater' section, Table 1):
% 1 or 'co2' = carbon dioxide CO2
% 2 or 'o2'  = oxygen O2
% Written by Matthew P. Humphreys, last updated 2016-11-06

switch lower(gas)
    case {1 'co2'}  % CO2, seawater salinity = 35
        A = 2116.8; B = -136.25; C = 4.7353; D = -0.092307; E = 0.0007555;
    case {2 'o2'} % O2, seawater salinity = 35
        A = 1920.4; B = -135.60; C = 5.2122; D = -0.10939; E = 0.00093777;
end %switch

sch = A + B*temp + C*temp.^2 + D*temp.^3 + E*temp.^4;
sch_valid = temp > -2 & temp < 40;
end %function co2f_schmidt14