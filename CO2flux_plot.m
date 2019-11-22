% Set inputs
u10 = (0:0.1:20)';
temp = 25;

% Calculate k for all references
refs = {'KRP06' 'KRP06-int' 'NML00' 'NML00-NS' 'RBHW19' 'TSW09' ...
    'W92' 'W14' 'W14-int'};
k = NaN(numel(u10), numel(refs));
for R = 1:numel(refs)
    k(:, R) = CO2flux_k_gasex(temp, u10, 'CO2', refs{R});
end % for R

% Draw figure
figure(1); clf; hold on
for R = 1:numel(refs)
    plot(u10, k(:, R))
end % for R
legend(refs, 'Location', 'nw')
xlabel('Wind speed / m/s')
ylabel('k / cm/hr')
