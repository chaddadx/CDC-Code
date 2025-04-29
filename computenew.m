%Plot where the NE is in pure strategies , depending on the values of alpha12 and alpha21.
%the Nash equilibrium can only be pure if there's a mutual best response
%check(min,min), (min,max), (max,min), (max,max)
%Define ranges for alpha12 and alpha21
% Initialize a grid to store results (e.g., 0 for pure).

% Parameters
T =7;          % Number of time steps
k = 3;          % (gamma_max = 1/k) for theorm 3
lambda =80;    % Cost coefficient

gamma_min = 0.01;
gamma_max =0.5;

% Define ranges for alpha12 and alpha21 
alpha12_range = linspace(0, 0.2, 1000);
alpha21_range = linspace(0, 0.2, 1000);


% store the NE type for each (alpha12, alpha21) pair.
% NE_type: 0 = No NE found, 1 = Pure NE
NE_type = zeros(length(alpha12_range), length(alpha21_range));
p_matrix = zeros(length(alpha12_range), length(alpha21_range));
q_matrix = zeros(length(alpha12_range), length(alpha21_range));


% Loop over alpha valuesc
for i = 1:length(alpha12_range)
    for j = 1:length(alpha21_range)
        alpha12 = alpha12_range(i);
        alpha21 = alpha21_range(j);
        
        % Compute utilities at the four strategy profiles
        U11 = compute_utility(gamma_min, gamma_min, T, k, lambda, alpha12, alpha21);
        U12 = compute_utility(gamma_min, gamma_max, T, k, lambda, alpha12, alpha21);
        U21 = compute_utility(gamma_max, gamma_min, T, k, lambda, alpha12, alpha21);
        U22 = compute_utility(gamma_max, gamma_max, T, k, lambda, alpha12, alpha21);
        
        % Utility matrices for each player
        % Rows correspond to Player 1's choices (1: gamma_min, 2: gamma_max)
        % Columns correspond to Player 2's choices (1: gamma_min, 2: gamma_max)
        U1 = [U11(1), U12(1);
            U21(1), U22(1)]; % Player 1's utilities
       
        
        U2 = [U11(2), U12(2);
            U21(2), U22(2)]; % Player 2's utilities
        
        % Check for pure Nash equilibrium
        pure_NE_found = false; %initially, it's set to false
        
        % Check each candidate pure strategy profile
        % Check candidate (1,1): (gamma_min, gamma_min)
        if U1(1,1) >= U1(2,1) && U2(1,1) >= U2(1,2)
            pure_NE_found = true;
        end
        
        
        % Check candidate (1,2): (gamma_min, gamma_max)
        if (U1(1,2) >= U1(2,2)) && (U2(1,2) >= U2(1,1))
            pure_NE_found = true;
        end
        
        
        % Check candidate (2,1): (gamma_max, gamma_min)
        if (U1(2,1) >= U1(1,1)) && (U2(2,1) >= U2(2,2))
            pure_NE_found = true;
        end
        
        
        % Check candidate (2,2): (gamma_max, gamma_max)
        if (U1(2,2) >= U1(1,2)) && (U2(2,2) >= U2(2,1))
            pure_NE_found = true; %If any of these conditions is true
        end
        
        
        
        if pure_NE_found
            NE_type(i,j) = 1;  % Mark as pure NE
        end
            
                NE_type(i,j) = 0;  % No valid NE detected
            end
        end


% PLOT THE RESULTS
figure;
imagesc(alpha21_range, alpha12_range, NE_type);
colorbar;
xlabel('\alpha_{21}');
ylabel('\alpha_{12}');
title('NE Type (0: No NE, 1: Pure NE, 2: Mixed NE)');
set(gca, 'YDir', 'normal');



%list all pure NEs for that alpha pair.
% Print the payoff matrices and NE information for each (alpha12, alpha21)
fprintf('Displaying payoff matrices and NE types for each (alpha12, alpha21) pair:\n');
for i = 1:length(alpha12_range)
    for j = 1:length(alpha21_range)
        alpha12 = alpha12_range(i);
        alpha21 = alpha21_range(j);
        
        
        U11 = compute_utility(gamma_min, gamma_min, T, k, lambda, alpha12, alpha21);
        U12 = compute_utility(gamma_min, gamma_max, T, k, lambda, alpha12, alpha21);
        U21 = compute_utility(gamma_max, gamma_min, T, k, lambda, alpha12, alpha21);
        U22 = compute_utility(gamma_max, gamma_max, T, k, lambda, alpha12, alpha21);
        
        U1 = [U11(1), U12(1);
            U21(1), U22(1)];
        
        
        U2 = [U11(2), U12(2);
            U21(2), U22(2)];
        
        fprintf('------------------------------\n');
        fprintf('alpha12 = %.3f, alpha21 = %.3f\n', alpha12, alpha21);
        fprintf('Player 1 Payoff Matrix:\n');
        fprintf('   [%.3f   %.3f]\n', U1(1,1), U1(1,2));
        fprintf('   [%.3f   %.3f]\n', U1(2,1), U1(2,2));
        fprintf('Player 2 Payoff Matrix:\n');
        fprintf('   [%.3f   %.3f]\n', U2(1,1), U2(1,2));
        fprintf('   [%.3f   %.3f]\n', U2(2,1), U2(2,2));
        
        % Displ%  ay NE
        if NE_type(i,j) == 1
            fprintf('NE Type: Pure NE\n');
            % Identify and print all pure NE candidate profiles
            NE_candidate = [];
            % Candidate (1,1): (gamma_min, gamma_min)
            if U1(1,1) >= U1(2,1) && U2(1,1) >= U2(1,2)
                NE_candidate = [NE_candidate; 1, 1];
            end
            
            % Candidate (1,2): (gamma_min, gamma_max)
            if U1(1,2) >= U1(2,2) && U2(1,2) >= U2(1,1)
                NE_candidate = [NE_candidate; 1, 2];
            end
            
            % Candidate (2,1): (gamma_max, gamma_min)
            if U1(2,1) >= U1(1,1) && U2(2,1) >= U2(2,2)
                NE_candidate = [NE_candidate; 2, 1];
            end
            
            % Candidate (2,2): (gamma_max, gamma_max)
            if U1(2,2) >= U1(1,2) && U2(2,2) >= U2(2,1)
                NE_candidate = [NE_candidate; 2, 2];
            end
            
            % Print each candidate in a readable format
            for c = 1:size(NE_candidate,1)
                cand = NE_candidate(c,:);
                if cand(1) == 1
                    gamma1_str = 'gamma_min';
                else
                    gamma1_str = 'gamma_max';
                end
                if cand(2) == 1
                    gamma2_str = 'gamma_min';
                else
                    gamma2_str = 'gamma_max';
                end
                fprintf('Pure NE candidate: (%s, %s) with payoffs: (%.3f, %.3f)\n', ...
                    gamma1_str, gamma2_str, U1(cand(1), cand(2)), U2(cand(1), cand(2)));
            end
        elseif NE_type(i,j) == 2
            fprintf('NE Type: Mixed NE, p = %.3f, q = %.3f\n', p_matrix(i,j), q_matrix(i,j));
        else
            fprintf('NE Type: No NE found\n');
        end
    end
end
fprintf('------------------------------\n');




% Create logical matrices for each candidate pure NE profile.
candidate11 = false(length(alpha12_range), length(alpha21_range)); % (gamma_min, gamma_min)
candidate12 = false(length(alpha12_range), length(alpha21_range)); % (gamma_min, gamma_max)
candidate21 = false(length(alpha12_range), length(alpha21_range)); % (gamma_max, gamma_min)
candidate22 = false(length(alpha12_range), length(alpha21_range)); % (gamma_max, gamma_max)

% Initialize NE_type to 0 
NE_type = zeros(length(alpha12_range), length(alpha21_range));

% Loop over the grid (using the same ranges as before)
for i = 1:length(alpha12_range)
    for j = 1:length(alpha21_range)
        alpha12 = alpha12_range(i);
        alpha21 = alpha21_range(j);
        
        % Compute utilities at the four strategy profiles
        U11 = compute_utility(gamma_min, gamma_min, T, k, lambda, alpha12, alpha21);
        U12 = compute_utility(gamma_min, gamma_max, T, k, lambda, alpha12, alpha21);
        U21 = compute_utility(gamma_max, gamma_min, T, k, lambda, alpha12, alpha21);
        U22 = compute_utility(gamma_max, gamma_max, T, k, lambda, alpha12, alpha21);
        
        % Form the payoff matrices (rows: Player 1, columns: Player 2)
        U1 = [U11(1), U12(1);
            U21(1), U22(1)];
        U2 = [U11(2), U12(2);
            U21(2), U22(2)];
        
        % Check candidate (1,1): (gamma_min, gamma_min)
        if U1(1,1) >= U1(2,1) && U2(1,1) >= U2(1,2)
            candidate11(i,j) = true;
            NE_type(i,j) = 1; % (min, min)
        end
        
        % Check candidate (1,2): (gamma_min, gamma_max)
        if U1(1,2) >= U1(2,2) && U2(1,2) >= U2(1,1)
            candidate12(i,j) = true;
            NE_type(i,j) = 2; % (min, max)
        end
        
        % Check candidate (2,1): (gamma_max, gamma_min)
        if U1(2,1) >= U1(1,1) && U2(2,1) >= U2(2,2)
            candidate21(i,j) = true;
            NE_type(i,j) = 3; % (max, min)
        end
        
        % Check candidate (2,2): (gamma_max, gamma_max)
        if U1(2,2) >= U1(1,2) && U2(2,2) >= U2(2,1)
            candidate22(i,j) = true;
            NE_type(i,j) = 4; % (max, max)
        end
    end
end



 colors = [...
    1.0  1.0  1.0;   % 1: White (min, min)
    0.6  0.8  1.0;   % 2: Light Blue (min, max)
    0.7  1.0  0.7;   % 3: Light Green (max, min)
    1.0  0.7  0.8];  % 4: Light Pink (max, max)


% Set the colormap
colormap(colors);




% Ensure NE_type contains values 1, 2, 3, and 4
% 1: (min, min), 2: (min, max), 3: (max, min), 4: (max, max)
imagesc(alpha21_range, alpha12_range, NE_type);

% Set the color axis to match the range of NE_type
caxis([1 4]);

% Add color bar with correct labels
%colorbar('Ticks', [1, 2, 3, 4], ... 'TickLabels', {'(\gamma_{min},\gamma_{min})', '(\gamma_{min},\gamma_{max})', '(\gamma_{max},\gamma_{min})', '(\gamma_{max},\gamma_{max})'});

% Label the axes and title
xlabel('\alpha_{21}');
ylabel('\alpha_{12}');


% Ensure the y-axis is in the correct direction
set(gca, 'YDir', 'normal');
% Save the figure as a PDF
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


% Suppose NE_type is the matrix of size length(alpha12_range) x length(alpha21_range)
% with values in {1,2,3,4} that you have computed above.

figure;
imagesc(alpha21_range, alpha12_range, NE_type);
set(gca, 'YDir', 'normal');
colormap(colors);  % your chosen colormap
caxis([1 4]);
xlabel('\alpha_{21}');
ylabel('\alpha_{12}');
hold on;

% Now pick 4 points interactively, one in each region:
disp('Click once inside each of the 4 NE regions ...');
[xClicks, yClicks] = ginput(4);

% Place text at those 4 clicked points:
text(xClicks(1), yClicks(1), '$(\gamma_{\min},\gamma_{\min})$', ...
    'Interpreter','latex','FontSize',12,'Color','k');
text(xClicks(2), yClicks(2), '$(\gamma_{\min},\gamma_{\max})$', ...
    'Interpreter','latex','FontSize',12,'Color','k');
text(xClicks(3), yClicks(3), '$(\gamma_{\max},\gamma_{\min})$', ...
    'Interpreter','latex','FontSize',12,'Color','k');
text(xClicks(4), yClicks(4), '$(\gamma_{\max},\gamma_{\max})$', ...
    'Interpreter','latex','FontSize',12,'Color','k');

hold off;

% Ensure the y-axis is in the correct direction
set(gca, 'YDir', 'normal');
% Save the figure as a PDF
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];



