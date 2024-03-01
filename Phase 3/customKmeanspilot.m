function [extractedPilotIndices] = customKmeanspilot(UE_X, UE_Y, K, pilotLength,R,L,sortedMaxIndices,C)
    % UE coordinates
    UE_coordinates = [UE_X, UE_Y];
    data = vertcat(UE_coordinates);
    
    % Number of UEs
    numUEs = K; 

    % Number of clusters
    numClusters = ceil(numUEs / pilotLength);
    k = numClusters;

    % Initial centroids
    centroids = datasample(data, numClusters, 'Replace', false);
    disp(centroids);
    
    % Initialize convergence flag
    converged = false;

    % Set convergence threshold
    convergenceThreshold = 0.001;

    % Initialize iteration counter
    iterations = 0;

    % Initialize cluster indices
    clusterIndices = zeros(numUEs, 1);

    % Initialize cell array to store cluster indices
    clusterUEIndices = cell(numClusters, 1);

%**************************************************************************
%****Cluster formation using K-means and handling 10 UEs in each cluster***
%**************************************************************************
    
    % Calculate distances from centroids
    distances = pdist2(data, centroids);
    
    % Assign points to clusters based on minimum distance
    [~, clusterIndices] = min(distances, [], 2);
    
    % Continue iterating until convergence
    while ~converged
        % Increment iteration count for performance check
        iterations = iterations + 1;
        
        % Update centroids
        oldCentroids = centroids;
        for i = 1:k
            clusterPoints = data(clusterIndices == i, :);
            centroids(i, :) = mean(clusterPoints);
        end
        
        % Calculate distances from centroids
        distances = pdist2(data, centroids);
        
        % Assign points to clusters based on minimum distance
        [~, clusterIndices] = min(distances, [], 2);
        
        % Check convergence by comparing centroid change
        centroidChange = sum(sqrt(sum((oldCentroids - centroids).^2, 2)));
        if centroidChange < convergenceThreshold
            converged = true;
        end
    end
    
    % Store cluster indices
    for i = 1:numClusters
        clusterUEIndices{i} = find(clusterIndices == i);
    end
    
    % Display cluster indices
    for i = 1:numel(clusterUEIndices)
        fprintf('ClusterIndices %d:\n', i);
        disp(clusterUEIndices{i});
    end
    disp("Number of iterations");
    disp(iterations);

    % Ensure each cluster holds at most 10 UEs
    for i = 1:numClusters
        clusterUEs = find(clusterIndices == i);
        centroid = centroids(i, :);
        disp(centroid);
        while numel(clusterUEs) > pilotLength
            % If more than 10 UEs, keep the 10 closest to the centroid
            distances = sqrt(sum((data(clusterUEs, :) - centroid).^2, 2));
            [~, sortedIdx] = sort(distances);
            clusterIndices(clusterUEs(sortedIdx(pilotLength+1:end))) = 0; % Set extra UEs to unassigned
            
            % Update clusterUEs after removing excess UEs
            clusterUEs = find(clusterIndices == i);
        end
    end

    % Reassign unassigned UEs to clusters
    unassignedUEs = find(clusterIndices == 0);
    for i = 1:numClusters
        numToAdd = min(pilotLength - nnz(clusterIndices == i), numel(unassignedUEs));
        if numToAdd > 0
            % Find closest unassigned UEs to the centroid
            centroid = centroids(i, :);
            distances = sqrt(sum((data(unassignedUEs, :) - centroid).^2, 2));
            [~, sortedIdx] = sort(distances);
            clusterIndices(unassignedUEs(sortedIdx(1:numToAdd))) = i; % Assign closest UEs
            
            % Update unassignedUEs after assigning UEs to clusters
            unassignedUEs = find(clusterIndices == 0);
        end
    end

    cluster_matrix = zeros(pilotLength, numClusters);
    % Display the final cluster indices
    for i = 1:numClusters
        clusterUEs = find(clusterIndices == i);
        cluster_matrix(:, i) = find(clusterIndices == i);
        fprintf('Indices in Cluster %d: %s\n', i, num2str(clusterUEs'));
    end

%**************************************************************************
%************************Pilot Assignment**********************************
%**************************************************************************
extractedPilotIndices = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C,cluster_matrix,1);
