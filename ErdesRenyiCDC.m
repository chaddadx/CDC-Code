% Parameters
N = 8000;              % Total number of nodes in the ER graph
T = 5;                 % Total time steps
num_iterations = 10000; % Number of simulation iterations for averaging



k = 5;
p = k / (N - 1);       % Connection probability for the ER graph

% Transition probabilities 
gamma1 = 0.4;
gamma2 = 0.35;
alpha12 = 0.15;
alpha21 = 0.2;

% Generate Erdős–Rényi random graph
% Create a symmetric adjacency matrix
A = rand(N) < p;       % random connections
A = triu(A, 1);        % keep upper triangular part
A = A + A';            % make symmetric
A(1:N+1:end) = 0;      % remove self-loops

% Simulation over multiple iterations
% State codes:
%    0: inactive (never activated)
%    1: active due to project 1
%    2: active due to project 2
%   -1: silent (already activated by project 1)
%   -2: silent (already activated by project 2)

average_cumulative_pc1 = zeros(1, T);
average_cumulative_pc2 = zeros(1, T); 

all_cumulative_pc1 = zeros(num_iterations, T);
all_cumulative_pc2 = zeros(num_iterations, T);

for iter = 1:num_iterations
    % Initialize the state matrix: rows = nodes, columns = time steps
    state = zeros(N, T);
    
    % At t=1, randomly choose 2 nodes (without replacement) and set:
    % one active with project 1 (state = 1) and one with project 2 (state = 2)
    nodes = randperm(N, 2);
    state(nodes(1), 1) = 1;
    state(nodes(2), 1) = 2;
    
    % Record cumulative activations at t=1
    cum_pc1 = zeros(1, T);
    cum_pc2 = zeros(1, T);
    cum_pc1(1) = sum(state(:,1) == 1 | state(:,1) == -1);
    cum_pc2(1) = sum(state(:,1) == 2 | state(:,1) == -2);
    
    % Time evolution for t = 2 to T
    for t = 2:T
        % Start with previous state
        new_state = state(:, t-1);
        
        % All nodes that were active at t-1 become silent at t.
        active_idx = find(state(:, t-1) == 1 | state(:, t-1) == 2);
        for i = 1:length(active_idx)
            idx = active_idx(i);
            if state(idx, t-1) == 1
                new_state(idx) = -1;
            elseif state(idx, t-1) == 2
                new_state(idx) = -2;
            end
        end
        
        % For each node that is still inactive (state 0), check if any active neighbor
        inactive_nodes = find(state(:, t-1) == 0);
        for j = 1:length(inactive_nodes)
            node = inactive_nodes(j);
            % Identify neighbors (indices where A(node,:) == 1)
            neighbors = find(A(node, :));
            % Consider only neighbors that were active (state 1 or 2) at time t-1
            active_neighbors = neighbors( state(neighbors, t-1) == 1 | state(neighbors, t-1) == 2 );
            % Randomize the order to remove any order bias
            if ~isempty(active_neighbors)
                active_neighbors = active_neighbors(randperm(length(active_neighbors)));
            end
            % Process each neighbor’s influence until the node is activated
            for n_idx = 1:length(active_neighbors)
                neighbor = active_neighbors(n_idx);
                if state(neighbor, t-1) == 1  % neighbor active with project 1
                    r = rand();
                    if r < gamma1
                        new_state(node) = 1;  % activated by project 1
                        break;
                    elseif r < gamma1 + alpha12
                        new_state(node) = 2;  % switched: activated by project 2
                        break;
                    end
                elseif state(neighbor, t-1) == 2  % neighbor active with project 2
                    r = rand();
                    if r < gamma2
                        new_state(node) = 2;  % activated by project 2
                        break;
                    elseif r < gamma2 + alpha21
                        new_state(node) = 1;  % switched: activated by project 1
                        break;
                    end
                end
            end
        end
        
        % Update state for time t
        state(:, t) = new_state;
        
        % Cumulative contributions: once a node is activated (active or silent) it remains so.
        cum_pc1(t) = sum(state(:, t) == 1 | state(:, t) == -1);
        cum_pc2(t) = sum(state(:, t) == 2 | state(:, t) == -2);
    end
    
    average_cumulative_pc1 = average_cumulative_pc1 + cum_pc1;
    average_cumulative_pc2 = average_cumulative_pc2 + cum_pc2;
    all_cumulative_pc1(iter, :) = cum_pc1;
    all_cumulative_pc2(iter, :) = cum_pc2;
end

% Compute averages over iterations
average_cumulative_pc1 = average_cumulative_pc1 / num_iterations;
average_cumulative_pc2 = average_cumulative_pc2 / num_iterations;

% Theoretical cumulative contributions
% Initial contributions

X1 = zeros(1,T);  % Theoretical cumulative contributions for Project 1
X2 = zeros(1,T);  % Theoretical cumulative contributions for Project 2
X1(1) = 1;
X2(1) = 1;
 


Z1 = zeros(1,T);
Z2 = zeros(1,T);
Z1(1) = 1;
Z2(1) = 1;

for t = 2:T

    Z1(t) = Z1(t-1) * k * gamma1 + Z2(t-1) * k * alpha21;
    Z2(t) =Z2(t-1) * k * gamma2 + Z1(t-1) * k * alpha12;
    X1(t) = X1(t-1) + Z1(t);
    X2(t) = X2(t-1) + Z2(t);
end

% Compute percentage errors between simulation and theory
pe_pc1 = mean( abs((X1 - average_cumulative_pc1) ./ X1) ) * 100;
pe_pc2 = mean( abs((X2 - average_cumulative_pc2) ./ X2) ) * 100;
disp(['Percentage Error (PE) for Project 1: ' num2str(pe_pc1) '%']);
disp(['Percentage Error (PE) for Project 2: ' num2str(pe_pc2) '%']);

% Visualization of the dynamics on the ER graph for one sample run

state_vis = zeros(N, T);
state_vis(:,1) = 0;
vis_nodes = randperm(N,2);
state_vis(vis_nodes(1),1) = 1;  % activate one node with project 1
state_vis(vis_nodes(2),1) = 2;  % activate one node with project 2

% Create graph object for plotting
G = graph(A);

figure;
for t = 1:T
    if t > 1
        new_state_vis = state_vis(:, t-1);
        % Active nodes become silent:
        active_idx = find(state_vis(:, t-1) == 1 | state_vis(:, t-1) == 2);
        for i = 1:length(active_idx)
            idx = active_idx(i);
            if state_vis(idx, t-1) == 1
                new_state_vis(idx) = -1;
            elseif state_vis(idx, t-1) == 2
                new_state_vis(idx) = -2;
            end
        end
        % Inactive nodes can be activated by active neighbors from previous time step.
        inactive_nodes = find(state_vis(:, t-1) == 0);
        for j = 1:length(inactive_nodes)
            node = inactive_nodes(j);
            neighbors = find(A(node, :));
            active_neighbors = neighbors( state_vis(neighbors, t-1) == 1 | state_vis(neighbors, t-1) == 2 );
            if ~isempty(active_neighbors)
                active_neighbors = active_neighbors(randperm(length(active_neighbors)));
            end
            for n_idx = 1:length(active_neighbors)
                neighbor = active_neighbors(n_idx);
                if state_vis(neighbor, t-1) == 1
                    r = rand();
                    if r < gamma1
                        new_state_vis(node) = 1;
                        break;
                    elseif r < gamma1 + alpha12
                        new_state_vis(node) = 2;
                        break;
                    end
                elseif state_vis(neighbor, t-1) == 2
                    r = rand();
                    if r < gamma2
                        new_state_vis(node) = 2;
                        break;
                    elseif r < gamma2 + alpha21
                        new_state_vis(node) = 1;
                        break;
                    end
                end
            end
        end
        state_vis(:, t) = new_state_vis;
    end
    
    % Plot the network at time t
    h = plot(G);
    node_colors = zeros(N, 3);
    for idx = 1:N
        s = state_vis(idx, t);
        if s == 1
            node_colors(idx, :) = [1 0 0];    % Red: active project 1
        elseif s == 2
            node_colors(idx, :) = [0 1 0];    % Green: active project 2
        elseif s == -1
            node_colors(idx, :) = [0 0 0];    % Black: silent (project 1)
        elseif s == -2
            node_colors(idx, :) = [1 1 0];    % Yellow: silent (project 2)
        else
            node_colors(idx, :) = [0 0 1];    % Blue: inactive
        end
    end
    h.NodeColor = node_colors;
    title(['ER Graph at time step t = ' num2str(t)]);
    pause(0.5);
end

% Plot simulated cumulative contributions vs. theoretical values
figure;
set(gcf, 'Position', [100, 100, 1000, 500]);

subplot(1, 2, 1);
plot(1:T, average_cumulative_pc1, '-r', 'LineWidth', 2, 'DisplayName', 'Simulated Z_1(t)');
hold on;
plot(1:T, X1, '--b', 'LineWidth', 2, 'DisplayName', 'Theoretical Z_1(t)');
xlabel('Time step t');
ylabel('Cumulative Contributions (Project 1)');
legend;

subplot(1, 2, 2);
plot(1:T, average_cumulative_pc2, '-k', 'LineWidth', 2, 'DisplayName', 'Simulated Z_2(t)');
hold on;
plot(1:T, X2, '--g', 'LineWidth', 2, 'DisplayName', 'Theoretical Z_2(t)');
xlabel('Time step t');
ylabel('Cumulative Contributions (Project 2)');
legend;



% Display results
disp(['PC1 Percentage Error: ', num2str(pe_pc1), '%']);
disp(['PC2 Percentage Error: ', num2str(pe_pc2), '%']);

disp('Final Cumulative Counts:');

disp(['Theoretical PC1: ', num2str(X1())]);
disp(['Simulated PC1: ', num2str(average_cumulative_pc1())]);


disp(['Theoretical PC2: ', num2str(X2())]);
disp(['Simulated PC2: ', num2str(average_cumulative_pc2())]);

