function [k,k_unc,k_valid] = co2f_k_gasex(temp,u10,gas,ref)
%co2f_k_gasex Calculates gas exchange coefficient <k> in cm/hr
% Part of co2flux() toolbox - https://github.com/mvdh7/co2flux
% Inputs: <temp> = seawater temperature / degC
%         <u10>  = wind speed at 10 m / m/s
%         <gas>  = name of gas ('co2')
%         <ref>  = reference to use for calculation ('w14')
% Written by Matthew P. Humphreys [last updated 2016-11-07]

switch lower(ref)
    
    case 't09'
        % Takahashi et al., 2009, Deep-Sea Res Pt II 56, 554-577
        [sch,sch_valid] = co2f_schmidt92(temp,gas);
        
        k = 0.26 * u10.^2 .* sqrt(660./sch); % cm/hr
        k_unc = 0.3 * k; % cm/hr
        
        k_valid = sch_valid;
    
    case 'w14' 
        % Wanninkhof, 2014, Limnol Oceanogr Methods 12, 351-362
        [sch,sch_valid] = co2f_schmidt14(temp,gas);

        k = 0.251 * u10.^2 .* sqrt(660./sch); % cm/hr
        k_unc = 0.2 * k; % cm/hr

        k_valid = sch_valid;
        
end %switch

end %function co2f_k_gasex