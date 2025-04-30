
% PARAMETERS
T = 5;
gamma_min = 0;
gamma_max = 0.3;
alpha12 = 0.1;
alpha21 = 0.1;

% Define ranges for lambda and k (adjust ranges as needed)
lambda_range = linspace(0, 50, 1000);
k_range = 2:8;  % Example range for k

% Preallocate matrix to store the candidate NE type for each (k, lambda) pair.
% Rows correspond to different k values and columns to lambda values.
NE_type = zeros(length(k_range), length(lambda_range));

% COMPUTE PURE NASH EQUILIBRIA OVER THE (k, lambda) GRID
for i = 1:length(k_range)
    k_val = k_range(i);
    for j = 1:length(lambda_range)
        lambda_val = lambda_range(j);
        
        % Compute utilities for the four strategy profiles:
        % (min,min), (min,max), (max,min), (max,max)
        U11 = compute_utility(gamma_min, gamma_min, T, k_val, lambda_val, alpha12, alpha21);
        U12 = compute_utility(gamma_min, gamma_max, T, k_val, lambda_val, alpha12, alpha21);
        U21 = compute_utility(gamma_max, gamma_min, T, k_val, lambda_val, alpha12, alpha21);
        U22 = compute_utility(gamma_max, gamma_max, T, k_val, lambda_val, alpha12, alpha21);
        
        % Build payoff matrices for Player 1 and Player 2
        U1 = [U11(1), U12(1); U21(1), U22(1)];
        U2 = [U11(2), U12(2); U21(2), U22(2)];
        
        % Determine pure NE candidate(s) only.
        candidates = [];
        if (U1(1,1) >= U1(2,1)) && (U2(1,1) >= U2(1,2))
            candidates(end+1) = 1;  % (min,min)
        end
        if (U1(1,2) >= U1(2,2)) && (U2(1,2) >= U2(1,1))
            candidates(end+1) = 2;  % (min,max)
        end
        if (U1(2,1) >= U1(1,1)) && (U2(2,1) >= U2(2,2))
            candidates(end+1) = 3;  % (max,min)
        end
        if (U1(2,2) >= U1(1,2)) && (U2(2,2) >= U2(2,1))
            candidates(end+1) = 4;  % (max,max)
        end
        
        % If no pure NE is found, mark as 0.
        if isempty(candidates)
            candidate_code = 0;
        else
            candidate_code = candidates(1); % use the first candidate found
        end
        
        NE_type(i, j) = candidate_code;
    end
end

% PLOT THE NASH EQUILIBRIA TYPES
figure;
% Define a colormap:
% Row 1: Candidate 0 (No pure NE) --> Gray.
% Rows 2-5: Pure NE candidates 1:(min,min), 2:(min,max), 3:(max,min), 4:(max,max)
colors = [0.8, 0.8, 0.8;   % 0: Gray (No pure NE)
          1.0, 1.0, 1.0;   % 1: White (min, min)
          0.6, 0.8, 1.0;   % 2: Light Blue (min, max)
          0.7, 1.0, 0.7;   % 3: Light Green (max, min)
          1.0, 0.7, 0.8];  % 4: Light Pink (max, max)

colormap(colors);
imagesc(lambda_range, k_range, NE_type);
set(gca, 'YDir', 'normal');
% Set color axis limits to cover codes 0 to 4.
caxis([-0.5 4.5]);
colorbar('Ticks', 0:4, 'TickLabels', {'None','(min,min)','(min,max)','(max,min)','(max,max)'});

xlabel('\mu');
ylabel('k');


% SAVE THE FIGURE AS A PDF
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, 'Pure_NE_Type_k_lambda.pdf', '-dpdf');


% For k = 3, find the corresponding row index in k_range
k_index = find(k_range == 3);

% Extract the NE type vector for k = 3
NE_k3 = NE_type(k_index, :);

% Find the first index where the Nash equilibrium type changes (indicating a threshold in mu)
change_idx = find(diff(NE_k3) ~= 0, 1, 'first');

% Extract the corresponding mu (lambda) value from lambda_range
mu_threshold = lambda_range(change_idx);

fprintf('For k = 3, the critical mu value is approximately %g.\n', mu_threshold);

