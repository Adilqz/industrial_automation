clear;
clc;
vrp_solver('XML100_1111_02.vrp');

function vrp_solver(file_name)
    % Parse VRP file
    [coords, demands, capacity] = parse_vrp_file(file_name);

    % Number of nodes
    num_nodes = size(coords, 1);

    % Calculate distance matrix
    dist_matrix = calc_distance_matrix(coords);

    % Solve VRP using a simple heuristic (Nearest Neighbor)
    [solution, route_costs] = solve_vrp_nn(dist_matrix, demands, capacity);

    % Calculate total cost
    total_cost = sum(route_costs);

    % Display the solution and costs in the command window
    disp('Solution:');
    disp(solution);
    disp('Route Costs:');
    for i = 1:length(route_costs)
        disp(['Route ', num2str(i), ': Cost = ', num2str(route_costs(i))]);
    end
    disp(['Total Cost: ', num2str(total_cost)]);

    % Plot the solution
    plot_solution(coords, solution);
end

function [coords, demands, capacity] = parse_vrp_file(file_name)
    fileID = fopen(file_name, 'r');
    tline = fgetl(fileID);
    
    while ischar(tline)
        if contains(tline, 'CAPACITY')
            capacity = sscanf(tline, 'CAPACITY : %d');
        elseif contains(tline, 'NODE_COORD_SECTION')
            break;
        end
        tline = fgetl(fileID);
    end

    coords = [];
    while ischar(tline) && ~contains(tline, 'DEMAND_SECTION')
        if ~contains(tline, 'NODE_COORD_SECTION')
            data = sscanf(tline, '%d %f %f');
            coords = [coords; data(2:3)'];
        end
        tline = fgetl(fileID);
    end
    assignin('base','coordinates',coords);

    demands = [];
    while ischar(tline) && ~contains(tline, 'DEPOT_SECTION')
        if ~contains(tline, 'DEMAND_SECTION')
            data = sscanf(tline, '%d %d');
            demands = [demands; data(2)];
        end
        tline = fgetl(fileID);
    end

    fclose(fileID);
end

function dist_matrix = calc_distance_matrix(coords)
    num_nodes = size(coords, 1);
    dist_matrix = zeros(num_nodes);

    for i = 1:num_nodes
        for j = 1:num_nodes
            dist_matrix(i, j) = sqrt((coords(i, 1) - coords(j, 1))^2 + (coords(i, 2) - coords(j, 2))^2);
        end
    end
    assignin('base','dist_matrix',dist_matrix);
end

function [solution, route_costs] = solve_vrp_nn(dist_matrix, demands, capacity)
    num_nodes = size(dist_matrix, 1);
    visited = false(1, num_nodes);
    current_node = 1;
    visited(1) = true;

    solution = [];
    current_route = [];
    current_load = 0;
    route_costs = [];
    current_route_cost = 0;

    while any(~visited)
        next_node = find_next_node(dist_matrix, current_node, visited, demands, current_load, capacity);

        if isempty(next_node)
            % Complete the route back to depot
            solution = [solution; current_route, 1];
            route_costs = [route_costs, current_route_cost + dist_matrix(current_node, 1)];
            current_route = [];
            current_load = 0;
            current_route_cost = 0;
            current_node = 1;
        else
            visited(next_node) = true;
            current_route = [current_route, next_node];
            current_load = current_load + demands(next_node);
            current_route_cost = current_route_cost + dist_matrix(current_node, next_node);
            current_node = next_node;
        end
    end

    if ~isempty(current_route)
        % Complete the last route back to depot
        solution = [solution; current_route, 1];
        route_costs = [route_costs, current_route_cost + dist_matrix(current_node, 1)];
    end
end

function next_node = find_next_node(dist_matrix, current_node, visited, demands, current_load, capacity)
    num_nodes = size(dist_matrix, 1);
    min_dist = inf;
    next_node = [];

    for i = 2:num_nodes  % Skip the depot (node 1)
        if ~visited(i) && (current_load + demands(i) <= capacity) && (dist_matrix(current_node, i) < min_dist)
            min_dist = dist_matrix(current_node, i);
            next_node = i;
        end
    end
end

function plot_solution(coords, solution)
    figure;
    hold on;
    grid on;

    % Plot all nodes
    plot(coords(:, 1), coords(:, 2), 'ko', 'MarkerFaceColor', 'k');
    text(coords(1, 1), coords(1, 2), 'Depot', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    for i = 2:size(coords, 1)
        text(coords(i, 1), coords(i, 2), num2str(i - 1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end

    % Plot routes
    colors = lines(size(solution, 1));
    for i = 1:size(solution, 1)
        route = solution(i, :);
        route_coords = coords(route, :);
        plot(route_coords(:, 1), route_coords(:, 2), '-', 'Color', colors(i, :), 'LineWidth', 2);
    end

    % Display total cost
    title('VRP Solution');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    hold off;
end



