% We simulates a branching process for two projects over T time
%A snapshot is a read copy of the system's data at a particular instant. 
% Saves the state of nodes at each time step for the final iteration .
% A growing “tree” that evolves over a fixed number of time steps.

% PARAMETERS
T =5;                   % Time steps 
num_iterations = 10000;   % Number of iterations
k = 5;             % Poisson parameter for number of children
gamma1 = 0.4;       
gamma2 = 0.35;      
alpha12 = 0.15;         % Probability to switch from Project 1 to Project 2 (child gets state 2)
alpha21 = 0.2;         % Probability to switch from Project 2 to Project 1 (child gets state 1)


% THEORETICAL EXPECTATIONS

% At t=1 each project starts with one active agent.

Z1 = zeros(1,T);  % Theoretical cumulative contributions for Project 1
Z2 = zeros(1,T);  % Theoretical cumulative contributions for Project 2
Z1(1) = 1;
Z2(1) = 1;
 

X1 = zeros(1,T);
X2 = zeros(1,T);
X1(1) = 1;
X2(1) = 1;

for t = 2:T

    X1(t) = X1(t-1) * k * gamma1 + X2(t-1) * k * alpha21;
    X2(t) =X2(t-1) * k * gamma2 + X1(t-1) * k * alpha12;
    Z1(t) = Z1(t-1) + X1(t);
    Z2(t) = Z2(t-1) + X2(t);
end


% SIMULATION OF THE BRANCHING PROCESS 
% Each node is stored as a row:
%   [global_id, parent_id, state, time, local_id, project]


% local_id: sequential number within the project .
% state:  1 (active project1), -1 (silent project1),
%         2 (active project2), -2 (silent project2),
%         0 (inactive).


% Inactive nodes will be colored blue. 
%Nodes are labeled as 'pc1-<local_id>' or 'pc2-<local_id>'.


% Edges are stored as rows: [parent_global_id, child_global_id].

%These will store the cumulative contributions from each iteration
%To compute an average at the end.
average_cum1 = zeros(1, T);
average_cum2 = zeros(1, T);


% A cell array to save the state (nodes) at each t for the final iteration. 
%This is used  for visualization.
lastIterationSnapshots = cell(T,1);  %A cell array where each cell will store the state of the network (the nodes) up to t for the final iteration.
% Parent–child relationships from the final iteration.
lastIterationEdges = []; %this will store the parent–child relationships (edges) from the final iteration.

for iter = 1:num_iterations
    %A matrix where each row represents an agent with six columns
    nodes = [];   % [global_id, parent_id, state, time, local_id, project]
    
    %Store pairs [parent_global_id, child_global_id].
    edges = [];   %  An empty matrix to record relationships [parent, child]
    next_id = 1;  % Global node counter to assign a unique global ID to each new node.
    
    % Track of the sequential numbering of agents within each project for the current iteration.
    %Count the nodes within each project so that each node gets a sequential “local” ID.
    local_counter1 = 0;   % For Project 1
    local_counter2 = 0;   % For Project 2
    
    %Create root nodes at t=1, two roots are created at t=1
    % Root for Project 1:
    %The node gets state 1 (active), local id is 1, and belongs to project 1.
    local_counter1 = local_counter1 + 1;
    nodes = [nodes; next_id, 0, 1, 1, local_counter1, 1]; %parent_id is set to 0, state is set to 1,Project number is 1
    root1 = next_id;
    next_id = next_id + 1;
    
    % Root for Project 2:
    local_counter2 = local_counter2 + 1;
    nodes = [nodes; next_id, 0, 2, 1, local_counter2, 2];
    root2 = next_id;
    next_id = next_id + 1;
    
    % A vector holding the global IDs of nodes that are currently “active”
    active = [root1, root2]; %the Active List for the First Time Step
    
    % Initialize cumulative contributions for this iteration:
    %Track the cumulative contributions at each t
    cum1 = zeros(1, T);
    cum2 = zeros(1, T);
    cum1(1) = 1;
    cum2(1) = 1;
    
    

    % For the last iteration, store at t = 1, for visualization
    if iter == num_iterations
        lastIterationSnapshots{1} = nodes(nodes(:,4)<=1, :); %for visualization purposes
    end
    
    %Simulation loop for time steps t = 2 to T 
    %the simulation updates the network by generating children nodes from the active nodes.
    for t = 2:T
        %A new list new_active is created to store the global IDs of nodes that will be active in the next time step.
        new_active = [];  
        for b = 1:length(active)
            parent_id = active(b);
            parent_row = nodes(nodes(:,1)==parent_id, :);
            parent_state = parent_row(3);
            parent_project = parent_row(6);
            
           
            
            %Generate Children for Each Parent  from a Poisson distribution with mean lambda
            nchild = poissrnd(k); 
            for j = 1:nchild
                %A random number r is generated.
                r = rand();
                if parent_state == 1  % (project 1)
                    if r < gamma1
                        child_state = 1;  % remains in Project 1
                    elseif r < gamma1 + alpha12
                        child_state = 2;  % switches to Project 2
                    else
                        child_state = 0;  % inactive 
                    end
                elseif parent_state == 2  % (project 2)
                    if r < gamma2
                        child_state = 2;  % remains in Project 2
                    elseif r < gamma2 + alpha21
                        child_state = 1;  % switches to Project 1
                    else
                        child_state = 0;  % inactive 
                    end
                end
                
                %  Determine the child's project and assign local id
                %If the child’s state is 1, it belongs to Project 1 in counting.
                if child_state == 1
                    project = 1;
                    local_counter1 = local_counter1 + 1;
                    local_id = local_counter1;
                    
                  %If the child’s state is 2, it belongs to Project 2.
                elseif child_state == 2
                    project = 2;
                    local_counter2 = local_counter2 + 1;
                    local_id = local_counter2;
                else  % inactive
                    project = parent_project;
                    if project == 1
                        local_counter1 = local_counter1 + 1;
                        local_id = local_counter1;
                    else
                        local_counter2 = local_counter2 + 1;
                        local_id = local_counter2;
                    end
                end
                
                % Create the child node: [global_id, parent_id, state, time, local_id, project]
                nodes = [nodes; next_id, parent_id, child_state, t, local_id, project];
                edges = [edges; parent_id, next_id];
                
                if child_state == 1 || child_state == 2
                    new_active(end+1) = next_id;
                end
                next_id = next_id + 1;
            end
            % Mark the parent as silent:
            %For a parent in Project 1, its state becomes -1.
            if parent_state == 1
                nodes(nodes(:,1)==parent_id, 3) = -1; %Column 3: State ,3 here refers to the 3 column of the nodes matrix, which is where the node's state is stored.
            %For a parent in Project 2, its state becomes -2.
            elseif parent_state == 2
                nodes(nodes(:,1)==parent_id, 3) = -2;
            end
        end
        
        %The active nodes list is updated to the new set.
        active = new_active;
        
        % Calculate contributions at time t :
        if ~isempty(new_active)
            %the third column of the nodes matrix.
            current_states = nodes(ismember(nodes(:,1), new_active), 3);
            X1_current = sum(current_states == 1);
            X2_current = sum(current_states == 2);
        else
            X1_current = 0;
            X2_current = 0;
        end
        cum1(t) = cum1(t-1) + X1_current;
        cum2(t) = cum2(t-1) + X2_current;
        
        % For the last iteration, store a snapshot of nodes up to time t:
        if iter == num_iterations
            lastIterationSnapshots{t} = nodes(nodes(:,4) <= t, :);
        end
    end
    
    average_cum1 = average_cum1 + cum1;
    average_cum2 = average_cum2 + cum2;
    
    % For the last iteration, store the complete edge list:
    if iter == num_iterations
        lastIterationEdges = edges; %Save the Edge List for the Final Iteration
    end
end

% Average cumulative contributions over all iterations:
average_cum1 = average_cum1 / num_iterations;
average_cum2 = average_cum2 / num_iterations;

% FINAL VISUALIZATION: Display network graphs for the last iteration from t = 1 to T
for t = 1:T
    figure;
    nodes_t = lastIterationSnapshots{t};
    
   
    valid_edges = [];  % include edges among nodes that have already been created.
    for i = 1:size(lastIterationEdges,1)
        child_id = lastIterationEdges(i,2);
        if any(nodes_t(:,1) == child_id)
            valid_edges = [valid_edges; lastIterationEdges(i,:)];
        end
    end
    if isempty(valid_edges)
        G = digraph; %directed graph or  it can be undirected using "graph"
        G = addnode(G, arrayfun(@num2str, nodes_t(:,1), 'UniformOutput', false));
    else
        G = digraph(valid_edges(:,1), valid_edges(:,2));
    end
    
    % Prepare labels: label each node as 'pc1-<local_id>' or 'pc2-<local_id>'
    labels = cell(size(nodes_t,1),1);
    for i = 1:size(nodes_t,1)
        proj = nodes_t(i,6);
        loc = nodes_t(i,5);
        if proj == 1
            labels{i} = sprintf('pc1-%d', loc);
        elseif proj == 2
            labels{i} = sprintf('pc2-%d', loc);
        end
    end
    
    % Determine node colors:
    % state  1  -> red,   state -1 -> black (Project 1)
    % state  2  -> green, state -2 -> yellow (Project 2)
    % state  0  -> blue (inactive)
    
    node_colors = zeros(size(nodes_t,1),3);
    for i = 1:size(nodes_t,1)
        st = nodes_t(i,3);
        if st == 1
            node_colors(i,:) = [1 0 0];   % red
        elseif st == -1
            node_colors(i,:) = [0 0 0];   % black
        elseif st == 2
            node_colors(i,:) = [0 1 0];   % green
        elseif st == -2
            node_colors(i,:) = [1 1 0];   % yellow
        elseif st == 0
            node_colors(i,:) = [0 0 1];   % blue for inactive
        end
    end
    
    h = plot(G, 'Layout', 'layered', 'Direction', 'right');
   %h = plot(G, 'Layout', 'layered'); % for the undirected
    title(sprintf('Network Graph at t = %d (Last Iteration)', t));
    h.NodeLabel = labels;
    h.NodeColor = node_colors;
    pause(1);
end

% Simulation vs. Theoretical 
figure;
set(gcf, 'Position', [100, 100, 1200, 500]);

subplot(1,2,1);
plot(1:T, average_cum1, '-r', 'LineWidth', 2, 'DisplayName', 'Simulated (Project 1)');
hold on;
plot(1:T, Z1, '--k', 'LineWidth', 2, 'DisplayName', 'Theoretical (Project 1)');
xlabel('Time step (t)');
ylabel('Cumulative Contributions');
title('Project 1 Contributions');
legend; grid on;

subplot(1,2,2);
plot(1:T, average_cum2, '-g', 'LineWidth', 2, 'DisplayName', 'Simulated (Project 2)');
hold on;
plot(1:T, Z2, '--k', 'LineWidth', 2, 'DisplayName', 'Theoretical (Project 2)');
xlabel('Time step (t)');
ylabel('Cumulative Contributions');
title('Project 2 Contributions');
legend; grid on;
% Save the figure as a PDF
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%Display the node table for inspection (for last iteration)
nodeTable = array2table(lastIterationSnapshots{T}, 'VariableNames', {'ID','Parent','State','Time','LocalID','Project'});
disp(nodeTable);

%Column 1: Global ID
%Column 2: Parent ID
%Column 3: State
%Column 4: Time
%Column 5: Local ID (agent number within its project)
%Column 6: Project (1 for Project 1, 2 for Project 2)

%percentage of error in comparison to the theoretical values
%the percentage error at the end
pe_pc1 = mean(abs((Z1 - average_cum1) ./ Z1)) * 100; % Percentage error for PC1
pe_pc2 = mean(abs((Z2 - average_cum2) ./Z2)) * 100; % Percentage error for PC2

% Display percentage errors
disp(['Percentage Error (PE) for PC1: ', num2str(pe_pc1), '%']);
disp(['Percentage Error (PE) for PC2: ', num2str(pe_pc2), '%']);

disp('Theoretical Cumulative Contributions : X̄ _1(t)'); disp(Z1);
disp('Simulated Cumulative Contributions (PC1): Z1'); disp(average_cum1);
disp('Theoretical Cumulative Contributions : X̄ _2(t)'); disp(Z2);
disp('Simulated Cumulative Contributions (PC2): Z2'); disp(average_cum2);
disp('Comparison: Theoretical vs. Simulated values are close when lambda approximates the average children per node.');





