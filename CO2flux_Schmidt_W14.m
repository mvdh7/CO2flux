function [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas)
%CO2flux_Schmidt_W14 Calculate Schmidt number following Wanninkhof (2014).
%
% Part of the CO2flux toolbox [https://github.com/mvdh7/co2flux].
% Written by Matthew P. Humphreys [last updated 2019-11-21].
%
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
    case {1 'co2'}
        % Carbon dioxide at practical salinity = 35.
        A = 2116.8;
        B = -136.25;
        C = 4.7353;
        D = -0.092307;
        E = 0.0007555;
    case {2 'o2'}
        % Oxygen at practical salinity = 35.
        A = 1920.4;
        B = -135.60;
        C = 5.2122;
        D = -0.10939;
        E = 0.00093777;
end
Sch = A + B*temp + C*temp.^2 + D*temp.^3 + E*temp.^4;
Sch_valid = temp > -2 & temp < 40;
