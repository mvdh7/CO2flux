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
%   ref = Reference to use for calculation (case insensitive):
%     'K06' = Krakauer et al. (2006), Tellus B 58(5), 390-471.
%     'K06i' = Krakauer et al. (2006) with a non-zero intercept, as
%              described by Ribas-Ribas et al. (2019).
%     'R19' = Ribas-Ribas et al. (2019), Geosciences 9(5), 230.
%     'T09' = Takahashi et al. (2009), Deep-Sea Res. Pt. II 56, 554-577.
%     'W14' = Wanninkhof (2014), Limnol. Oceanogr. Methods 12, 351-362.
%     'W14i' = Wanninkhof (2014) with a non-zero intercept, as described by
%              Ribas-Ribas et al. (2019).
% Outputs:
%   k = Gas exchange coefficient in cm/hr.
%   k_unc = Uncertainty in k in cm/hr.
%   k_valid = Is the input temperature within the valid range?
switch lower(ref)
    case 't09'
        % Takahashi et al. (2009), Deep-Sea Res. Pt. II 56, 554-577.
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = 0.26*u10.^2.*sqrt(660./Sch); % cm/hr
        k_unc = 0.3*k; % cm/hr
        k_valid = Sch_valid;
    case 'w14' 
        % Wanninkhof (2014), Limnol. Oceanogr. Methods 12, 351-362.
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = 0.251*u10.^2.*sqrt(660./Sch); % cm/hr
        k_unc = 0.2*k; % cm/hr
        k_valid = Sch_valid;
    case 'k06'
        % Krakauer et al. (2006), Tellus B 58(5), 390-471.
        % via Ribas-Ribas et al. (2019)
        [Sch, Sch_valid] = CO2flux_Schmidt_W92(temp, gas);
        k = 2.275*u10.*sqrt(660./Sch); % cm/hr
        k_unc = NaN*k;
        k_valid = Sch_valid;
    case {'w14i' 'mr2'}
        % Ribas-Ribas et al. (2019) "Wanninkhof + intercept"
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = (10.7 + 0.30*u10.^2).*sqrt(660./Sch); % cm/hr
        k_unc = NaN*k;
        k_valid = Sch_valid;
    case {'r19' 'msf'}
        % Mustaffa et al., in prep
        % Equation is pers. comm. from Mariana Ribas-Ribas
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = (10.71 + 0.06*u10.^2).*sqrt(660./Sch); % cm/hr
        k_unc = NaN * k;
        k_valid = Sch_valid;
    case {'k06i' 'bern'}
        % Ribas-Ribas et al. (2019) "Krakauer + intercept"
        [Sch, Sch_valid] = CO2flux_Schmidt_W14(temp, gas);
        k = (11 + 2.275*u10).*sqrt(660./Sch); % cm/hr
        k_unc = NaN*k;
        k_valid = Sch_valid;
end
