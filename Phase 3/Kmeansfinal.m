function [extractedPilotIndices] = Kmeansfinal(UE_X, UE_Y, K, pilotLength,R,L,sortedMaxIndices,C)
    
    % UE coordinates
    UE_coordinates = [UE_X, UE_Y];
    data = vertcat(UE_coordinates);
    
    numUEs = K;
    % Number of clusters
    numClusters = ceil(numUEs / pilotLength);
    k = numClusters;

    % Initial centroids
    centroids = datasample(data, numClusters, 'Replace', false);
    %disp(centroids);
    
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

% Store converged centroids
convergedCentroids = centroids;

% Calculate distances of UEs from converged centroids
UE_Distances = pdist2(data, convergedCentroids);

disp("convergedcentroids")
disp(convergedCentroids);
disp("UE Distances");
disp(UE_Distances);

% Clear clusterIndices for reassignment based on distances
clusterIndices = zeros(numUEs, 1);

% Assign UEs based on distances from centroids and check capacity
for iter = 1:numUEs
    minDist = Inf;
    minCluster = 0;
    minUE = 0;

    % Iterate through each UE
    for UE = 1:numUEs
        if clusterIndices(UE) == 0
            % Iterate through each cluster's centroid
            for cluster = 1:numClusters
                dist = UE_Distances(UE, cluster);
                if nnz(clusterIndices == cluster) < pilotLength && dist < minDist
                    minDist = dist;
                    minCluster = cluster;
                    minUE = UE;
                end
            end
        end
    end

    % Assign the UE to the cluster with minimum distance and available space
    if minCluster > 0
        clusterIndices(minUE) = minCluster;
    end
end

cluster_matrix = zeros(pilotLength, numClusters);
% Display cluster indices
for i = 1:numClusters
    clusterUEIndices = find(clusterIndices == i);
    cluster_matrix(:, i) = find(clusterIndices == i);
    fprintf('Indices in Cluster %d: %s\n', i, num2str(clusterUEIndices'));
end

%**************************************************************************
%************************Pilot Assignment**********************************
%**************************************************************************
%extractedPilotIndices = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C,cluster_matrix,1);

%% Random pilot assignment in each cluster 
extractedPilotIndices = zeros(numUEs, 1); 
pilotMatrix = zeros(pilotLength,numClusters);
for j = 1:numClusters
    %generate pilotLength random integers in the range of 1 to pilotLength, without reusing them.
    pilotMatrix(:, j) = randperm(pilotLength,pilotLength); 

    %store the pilot indices in the corresponding UE indices
    extractedPilotIndices(cluster_matrix(:,j)) = pilotMatrix(:,j); 
end

%cluster_matrix [p, K/p] stores which UE is in which cluster -
%cluster_matrix (5, 1) displays the 5th user that belongs to cluster 1

%extractedPilotIndices is essentially the pilotIndex, storing which UE has
%which pilot


extractedPilotIndices = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C,cluster_matrix,1,extractedPilotIndices);
