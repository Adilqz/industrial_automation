%% Main
clear;
filepath = 'C:\Users\adil2\Documents\MATLAB\Project 3\XML100_1111_01.vrp';
vrpData = readVRPFile(filepath);
plotVRPNodes(vrpData);
distanceMatrix = calculateDistanceMatrix(vrpData.nodeCoords);
[routes, totalDistance] = solveVRP(vrpData, distanceMatrix);
plotVRPSolution(vrpData, routes);

disp(['Total Distance: ', num2str(totalDistance)]);

%% Function for reading the file
function vrpData = readVRPFile(filepath)
    % Open the file
    fileID = fopen(filepath, 'r');
    if fileID == -1
        error('Cannot open file: %s', filepath);
    end

    % Initialize variables
    nodeCoords = [];
    demands = [];
    timeWindows = [];
    capacity = 0;
    numVehicles = 0;

    % Read file line by line
    while ~feof(fileID)
        line = fgetl(fileID);

        % Skip empty lines
        if isempty(line)
            continue;
        end

        % Parse the header information
        if startsWith(line, 'DIMENSION')
            dimension = sscanf(line, 'DIMENSION : %d');
        elseif startsWith(line, 'VEHICLES')
            numVehicles = sscanf(line, 'VEHICLES : %d');
        elseif startsWith(line, 'CAPACITY')
            capacity = sscanf(line, 'CAPACITY : %d');
        elseif startsWith(line, 'NODE_COORD_SECTION')
            nodeCoords = zeros(dimension, 3);
            for i = 1:dimension
                line = fgetl(fileID);
                nodeCoords(i, :) = sscanf(line, '%d %f %f');
            end
        elseif startsWith(line, 'DEMAND_SECTION')
            demands = zeros(dimension, 2);
            for i = 1:dimension
                line = fgetl(fileID);
                demands(i, :) = sscanf(line, '%d %d');
            end
        elseif startsWith(line, 'TIME_WINDOW_SECTION')
            timeWindows = zeros(dimension, 3);
            for i = 1:dimension
                line = fgetl(fileID);
                timeWindows(i, :) = sscanf(line, '%d %d %d');
            end
        elseif startsWith(line, 'DEPOT_SECTION')
            line = fgetl(fileID); % Read depot line
            depot = sscanf(line, '%d');
            % Read the end of depot section
            line = fgetl(fileID);
            if ~strcmp(line, '-1')
                error('Depot section format error');
            end
        end
    end

    % Close the file
    fclose(fileID);

    % Create output structure
    vrpData = struct();
    vrpData.nodeCoords = nodeCoords;
    vrpData.demands = demands;
    vrpData.timeWindows = timeWindows;
    vrpData.capacity = capacity;
    vrpData.numVehicles = numVehicles;
    vrpData.depot = depot;
end

%% Function for plotting the map
function plotVRPNodes(vrpData)
    % Extract the node coordinates and depot
    nodeCoords = vrpData.nodeCoords;
    depotIndex = vrpData.depot;
    
    % Separate depot and customer coordinates
    depotCoord = nodeCoords(depotIndex, 2:3);
    customerCoords = nodeCoords(:, 2:3);
    customerCoords(depotIndex, :) = []; % Remove depot coordinates from customers
    
    % Create the plot
    figure;
    hold on;
    
    % Plot customers
    plot(customerCoords(:, 1), customerCoords(:, 2), 'bo', 'MarkerSize', 5, 'DisplayName', 'Customers');
    
    % Plot depot
    plot(depotCoord(1), depotCoord(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Depot');
    
    % Add labels and legend
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('VRP Nodes and Depot');
    legend('show');
    
    hold off;
end

%% Calculation of distance matrix
function distanceMatrix = calculateDistanceMatrix(nodeCoords)
    numNodes = size(nodeCoords, 1);
    distanceMatrix = zeros(numNodes, numNodes);
    
    for i = 1:numNodes
        for j = 1:numNodes
            if i ~= j
                distanceMatrix(i, j) = sqrt((nodeCoords(i, 2) - nodeCoords(j, 2))^2 + ...
                                            (nodeCoords(i, 3) - nodeCoords(j, 3))^2);
            end
        end
    end
end

%% Solution using Nearest Neighbour Heuristic
function [routes, totalDistance] = solveVRP(vrpData, distanceMatrix)
    numNodes = size(vrpData.nodeCoords, 1);
    numVehicles = vrpData.numVehicles;
    capacity = vrpData.capacity;
    demands = vrpData.demands(:, 2);
    depot = vrpData.depot;

    remainingDemands = demands;
    routes = cell(numVehicles, 1);
    totalDistance = 0;

    for v = 1:numVehicles
        currentLoad = 0;
        currentRoute = depot;
        currentNode = depot;

        while true
            nextNode = findNextNode(currentNode, remainingDemands, currentLoad, capacity, distanceMatrix);
            if nextNode == depot || isempty(nextNode)
                break;
            end
            currentRoute = [currentRoute, nextNode];
            currentLoad = currentLoad + remainingDemands(nextNode);
            remainingDemands(nextNode) = 0;
            totalDistance = totalDistance + distanceMatrix(currentNode, nextNode);
            currentNode = nextNode;
        end
        
        currentRoute = [currentRoute, depot];
        totalDistance = totalDistance + distanceMatrix(currentNode, depot);
        routes{v} = currentRoute;
    end
end

function nextNode = findNextNode(currentNode, remainingDemands, currentLoad, capacity, distanceMatrix)
    numNodes = length(remainingDemands);
    minDistance = inf;
    nextNode = [];

    for i = 1:numNodes
        if i ~= currentNode && remainingDemands(i) > 0 && (currentLoad + remainingDemands(i)) <= capacity
            if distanceMatrix(currentNode, i) < minDistance
                minDistance = distanceMatrix(currentNode, i);
                nextNode = i;
            end
        end
    end

    if isempty(nextNode)
        nextNode = currentNode; % Return to depot if no next node is found
    end
end

%% Plot of the solution
function plotVRPSolution(vrpData, routes)
    nodeCoords = vrpData.nodeCoords;
    depotIndex = vrpData.depot;
    
    depotCoord = nodeCoords(depotIndex, 2:3);
    customerCoords = nodeCoords(:, 2:3);
    customerCoords(depotIndex, :) = [];

    figure;
    hold on;

    % Plot customers
    plot(customerCoords(:, 1), customerCoords(:, 2), 'bo', 'MarkerSize', 5, 'DisplayName', 'Customers');

    % Plot depot
    plot(depotCoord(1), depotCoord(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Depot');
    
    % Plot routes
    colors = lines(length(routes));
    for v = 1:length(routes)
        route = routes{v};
        routeCoords = nodeCoords(route, 2:3);
        plot(routeCoords(:, 1), routeCoords(:, 2), '-o', 'Color', colors(v, :), 'DisplayName', ['Vehicle ' num2str(v)]);
    end
    
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('VRP Solution');
    legend('show');

    hold off;
end


