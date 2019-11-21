function [FCO2, FCO2_unc, FCO2_valid, delpCO2, k, k0] = ...
    CO2flux(temp, psal, u10, pCO2air, pCO2air_unc, pCO2sw, pCO2sw_unc, ref)
%CO2flux Calculate air-sea CO2 flux with uncertainty propagation.
%
% Part of the CO2flux toolbox [https://github.com/mvdh7/co2flux].
% Written by Matthew P. Humphreys [last updated 2019-11-21].
%
% Inputs:
%   temp = Sea surface temperature, in degrees-C.
%   psal = Practical salinity.
%   u10 = Wind speed at 10 m above the sea surface in m/s.
%   pCO2air = Atmospheric pCO2 in micro-atm.
%   pCO2air_unc = Uncertainty in pCO2air in micro-atm.
%   pCO2sw = surface seawater pCO2 in micro-atm.
%   pCO2sw_unc = uncertainty in pCO2sw in micro-atm.
%   ref = choice of gas transfer coefficient equation:
%     'w14' = Wanninkhof (2014), Limnol. Oceanogr. Methods 12, 351-362.
%     't09' = Takahashi et al. (2009), Deep-Sea Res. Pt. II 56, 554-577.
%
% Outputs:
%   FCO2 = sea-to-air CO2 flux in micro-mol/(m^2*hr).
%   FCO2_unc = uncertainty in FCO2 in micro-mol/(m^2*hr).
%   FCO2_valid = are inputs temp & psal valid for Schmidt no. calculation?
[k, k_unc, k_valid] = CO2flux_k_gasex(temp, u10, 'co2', ref); % cm/hr
[k0, k0_unc, k0_valid] = CO2flux_HenrysCO2(temp, psal); % mol/(l*atm)
Tr = k.*k0*10; % mol/(m^2*hr*atm)
Tr_unc = sqrt((k_unc./k).^2 + (k0_unc./k0).^2).*Tr;
Tr_valid = k_valid & k0_valid;
delpCO2 = pCO2sw - pCO2air;
delpCO2_unc = sqrt(pCO2sw_unc.^2 + pCO2air_unc.^2);
FCO2 = delpCO2 .* Tr; % umol/(m^2*hr)
FCO2_unc = sqrt((delpCO2_unc./delpCO2).^2 + (Tr_unc./Tr).^2).*FCO2;
  % umol/(m^2*hr)
FCO2_valid = Tr_valid;
