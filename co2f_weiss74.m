function [k0,k0_c90,k0_weiss,k0_weiss_unc,k0_valid] = co2f_weiss74(t,s)
%co2f_weiss74 Calculates Henry's constant <k0> in mol/(l atm)
% Part of co2flux() toolbox - https://github.com/mvdh7/co2flux
% following Weiss (1974) Mar Chem 2, 203-215
% Inputs: <t> = temperature / degC
%         <s> = practical salinity

% Data from Table II of Weiss (1974)
tmp = [-1;0;1;2;3;4;5;6;8;10;12;14;16;18;20;22;24;26;28;30;32;34;36;38;40];
sal = [0 10 20 30 34 35 36 38 40];
[salg,tmpg] = meshgrid(sal,tmp);
k0_vol = [  NaN   NaN 7.273 6.903 6.760 6.724 6.689 6.620 6.551;
          7.758 7.364 6.990 6.635 6.498 6.465 6.431 6.364 6.298;
          7.458 7.081 6.723 6.382 6.251 6.219 6.187 6.123 6.060;
          7.174 6.813 6.469 6.143 6.017 5.986 5.955 5.894 5.833;
          6.905 6.558 6.229 5.916 5.795 5.766 5.736 5.677 5.619;
          6.650 6.317 6.001 5.701 5.585 5.557 5.528 5.472 5.416;
          6.408 6.088 5.785 5.497 5.386 5.358 5.331 5.277 5.223;
          6.178 5.871 5.580 5.303 5.196 5.170 5.144 5.092 5.040;
          5.751 5.469 5.200 4.945 4.846 4.822 4.797 4.749 4.702;
          5.366 5.105 4.857 4.621 4.529 4.507 4.485 4.440 4.396;
          5.017 4.776 4.546 4.327 4.243 4.222 4.201 4.160 4.119;
          4.700 4.477 4.264 4.062 3.983 3.964 3.945 3.906 3.869;
          4.412 4.205 4.008 3.820 3.747 3.729 3.712 3.676 3.641;
          4.149 3.958 3.775 3.600 3.533 3.516 3.499 3.466 3.434;
          3.910 3.732 3.562 3.400 3.337 3.322 3.306 3.275 3.245;
          3.691 3.526 3.368 3.217 3.158 3.144 3.130 3.101 3.073;
          3.491 3.337 3.190 3.050 2.995 2.982 2.968 2.942 2.915;
          3.307 3.164 3.027 2.897 2.846 2.833 2.821 2.796 2.771;
          3.138 3.005 2.878 2.756 2.709 2.697 2.685 2.662 2.639;
          2.983 2.859 2.741 2.627 2.583 2.572 2.561 2.540 2.518;
          2.840 2.725 2.615 2.509 2.468 2.457 2.447 2.427 2.407;
          2.708 2.601 2.498 2.400 2.361 2.352 2.342 2.323 2.305;
          2.587 2.487 2.391 2.299 2.263 2.254 2.246 2.228 2.211;
          2.474 2.382 2.292 2.207 2.173 2.165 2.157 2.140 2.124;
          2.370 2.284 2.201 2.121 2.090 2.082 2.074 2.059 2.044] * 1e-2;

% Generate fit
tkg = tmpg + 273.15;

k0_fit = fittype(['A1 + A2*(100/tk) + A3*log(tk/100)' ...
    '+ sal * (B1 + B2*(tk/100) + B3*(tk/100)^2)'], ...
    'independent',{'tk' 'sal'}, ...
    'coefficients',{'A1' 'A2' 'A3' 'B1' 'B2' 'B3'});
L = ~isnan(k0_vol);
k0_opts = fitoptions('method','nonlinearleastsquares');
k0_opts.Display = 'Off';
k0_opts.StartPoint = [-58 90 22 0 0 0];
k0_fit = fit([tkg(L) salg(L)],log(k0_vol(L)),k0_fit,k0_opts);

% Predict k0
tk = t + 273.15;
k0 = NaN(size(t));
k0(:) = exp(k0_fit([tk(:) s(:)]));
k0_c90 = NaN(size(t));
k0_c90(:) = diff(exp(predint(k0_fit,[tk(:) s(:)],0.9)),[],2);

% Calculate using Weiss (1974) coefficients
A1 = -58.0931;
A2 =  90.5069;
A3 =  22.2940;
B1 =   0.027766;
B2 =  -0.025888;
B3 =   0.0050578;
k0_weiss = exp(A1 + A2*(100./tk) + A3*log(tk/100) ...
    + s .* (B1 + B2*(tk/100) + B3*(tk/100).^2));

k0_weiss_unc = 1e-4 * ones(size(t));

% Validity logical
k0_valid = t >= -1 & t <= 40 & s >= 0 & s <= 40;

end %function co2f_weiss74