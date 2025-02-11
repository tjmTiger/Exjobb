function G = generate_directed_connected_graph_try(numNodes, numEdges)
    % Check for valid input
    if numEdges < numNodes
        error('The number of edges must be at least equal to the number of nodes to ensure connectivity.');
    end

    % Initialize directed edges
    edges = [];

    % Create a connected base graph using a directed chain
    for i = 1:numNodes-1
        edges = [edges; i, i+1];  % Create a directed edge from node i to i+1
    end

    % Randomly add remaining edges to meet the numEdges requirement
    while size(edges, 1) < numEdges
        % Generate random edge
        fromNode = randi(numNodes);
        toNode = randi(numNodes);
        
        % Ensure no self-loops and no duplicate edges
        if fromNode ~= toNode && ~ismember([fromNode, toNode], edges, 'rows')
            edges = [edges; fromNode, toNode];
        end
    end

    % Create the directed graph
    G = digraph(edges(:,1), edges(:,2));
    
    % Plot the directed graph
    figure;
    plot(G, 'Layout', 'layered');
    title('Directed Connected Graph');
end