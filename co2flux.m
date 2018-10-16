function [FCO2,FCO2_unc,FCO2_valid,delpCO2,k,k0_weiss] ...
    = co2flux(temp,sal,u10,pco2air,pco2air_unc,pco2sw,pco2sw_unc,ref)
%co2flux Calculates air-sea CO2 flux
% Part of co2flux() toolbox - https://github.com/mvdh7/co2flux
% OUTPUTS:
% -        FCO2: sea-to-air CO2 flux, in umol/(m^2 hr)
% -    FCO2_unc: uncertainty in <FCO2>, in umol/(m^2 hr)
% -  FCO2_valid: input <temp> & <sal> valid for Schmidt no. calculation?
% INPUTS:
% -        temp: sea surface temperature, in degrees C
% -         sal: practical salinity
% -         u10: wind speed at 10 m, in m/s
% -     pco2air: atmospheric pCO2, in uatm
% - pco2air_unc: uncertainty in <pco2air>, in uatm
% -      pco2sw: surface seawater pCO2, in uatm
% -  pco2sw_unc: uncertainty in <pco2sw>, in uatm
% -         ref: choice of gas transfer coefficient equation -
%        'w14' Wanninkhof, 2014, Limnol Oceanogr Methods 12, 351-362
%        't09' Takahashi et al., 2009, Deep-Sea Res Pt II 56, 554-577
% Written by Matthew P. Humphreys [last updated 2018-10-16]

[k,k_unc,k_valid] = co2f_k_gasex(temp,u10,'co2',ref); % cm/hr
[~,~,k0_weiss,k0_weiss_unc,k0_valid] = co2f_weiss74(temp,sal); % mol/(l atm)

Tr = k .* k0_weiss * 1e1; % mol / (m^2 hr atm)
Tr_unc = sqrt((k_unc./k).^2 + (k0_weiss_unc./k0_weiss).^2) .* Tr;
Tr_valid = k_valid & k0_valid;

delpCO2 = pco2sw - pco2air;
delpCO2_unc = sqrt(pco2sw_unc.^2 + pco2air_unc.^2);

FCO2 = delpCO2 .* Tr; % umol / (m^2 hr)
FCO2_unc = sqrt((delpCO2_unc./delpCO2).^2 + (Tr_unc./Tr).^2) .* FCO2;
  % umol / (m^2 hr)
FCO2_valid = Tr_valid;

end %function co2flux