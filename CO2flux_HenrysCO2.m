function [k0, k0_unc, k0_valid] = CO2flux_HenrysCO2(temp, psal)
%CO2flux_HenrysCO2 Calculate Henry's constant for CO2 in mol/(l*atm).
% Following Weiss (1974), Mar. Chem. 2, 203-215.
% Part of the CO2flux toolbox [https://github.com/mvdh7/CO2flux].
% Written by Matthew P. Humphreys [v1.0.0, last updated 2019-11-22].
% Inputs:
%   temp = Seawater temperature in degrees-C.
%   psal = Practical salinity.
% Outputs:
%   k0 = Henry's constant.
%   k0_unc = Uncertainty in k0.
%   k0_valid = Are the input temperature and salinity within their valid
%              ranges?
A1 = -58.0931;
A2 = 90.5069;
A3 = 22.2940;
B1 = 0.027766;
B2 = -0.025888;
B3 = 0.0050578;
tempK = temp + 273.15;
k0 = exp(A1 + A2*(100./tempK) + A3*log(tempK/100) + ...
    psal.*(B1 + B2*(tempK/100) + B3*(tempK/100).^2));
k0_unc = 1e-4*ones(size(temp));
k0_valid = temp >= -1 & temp <= 40 & psal >= 0 & psal <= 40;
