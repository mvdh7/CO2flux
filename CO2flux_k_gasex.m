function [k, k_unc, k_valid] = CO2flux_k_gasex(temp, u10, gas, ref)
%CO2flux_k_gasex Calculate gas exchange coefficient in cm/hr.
% Part of the CO2flux toolbox [https://github.com/mvdh7/CO2flux].
% Written by Matthew P. Humphreys [v1.0.0, last updated 2019-11-22].
% Inputs:
%   temp = Seawater temperature in degrees-C.
%   u10 = Wind speed at 10 m above the sea surface in m/s.
%   gas = Name of gas for Schmidt number calculation (case insensitive):
%     'CO2' = Carbon dioxide.
%     'O2' = Oxygen.
%   ref = Reference code to use for calculation (case insensitive).
%         Available options are (see CO2flux_citations function and
%         explanatory notes in the code below):
%     'KRP06', 'KRP06-int', 'NML00', 'NML00-NS', 'RBHW19', 'TSW09', 'W92',
%     'W14' 'W14-int'.
% Outputs:
%   k = Gas exchange coefficient in cm/hr.
%   k_unc = Uncertainty in k in cm/hr.
%   k_valid = Is the input temperature within the valid range?
% 
% TODO: allow user to select which Schmidt equation to use
% 
switch lower(ref)
    case {'tsw09' 't09'}
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = 0.26*u10.^2.*sqrt(660./Sch); % cm/hr
        k_unc = 0.3*k; % cm/hr
        k_valid = Sch_valid;
    case 'w14' 
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = 0.251*u10.^2.*sqrt(660./Sch); % cm/hr
        k_unc = 0.2*k; % cm/hr
        k_valid = Sch_valid;
    case 'w92'
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = 0.31*u10.^2.*sqrt(660./Sch); % cm/hr
        k_unc = 0.25*k; % cm/hr
        k_valid = Sch_valid;
    case {'krp06' 'k06'}
        % KRP06 via RBHW19.
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = 2.275*u10.*sqrt(660./Sch); % cm/hr
        k_unc = NaN*k;
        k_valid = Sch_valid;
    case {'w14-int' 'mr2'}
        % "Wanninkhof + intercept" equation of RHBW19.
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = (10.7 + 0.30*u10.^2).*sqrt(660./Sch); % cm/hr
        k_unc = NaN*k;
        k_valid = Sch_valid;
    case {'rbhw19' 'msf'}
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = (10.71 + 0.06*u10.^2).*sqrt(660./Sch); % cm/hr
        k_unc = NaN * k;
        k_valid = Sch_valid;
    case {'krp06-int' 'bern'}
        % "Krakauer + intercept" equation of RHBW19.
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = (11 + 2.275*u10).*sqrt(660./Sch); % cm/hr
        k_unc = NaN*k;
        k_valid = Sch_valid;
    case 'nml00-ns'
        % North-Sea-only version.
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = (0.23*u10.^2 + 0.1*u10).*sqrt(600./Sch);
        k_unc = 0.19*k;
        k_valid = Sch_valid;
    case 'nml00'
        % All-data version.
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = (0.222*u10.^2 + 0.333*u10).*sqrt(600./Sch);
        k_unc = 0.2*k;
        k_valid = Sch_valid;
end % switch
