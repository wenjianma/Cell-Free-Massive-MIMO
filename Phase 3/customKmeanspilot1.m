function [extractedPilotIndices] = customKmeanspilot1(UE_X, UE_Y, K, pilotLength,R,L,sortedMaxIndices,C)

    % UE coordinates
    UE_coordinates = [UE_X, UE_Y];
    data = vertcat(UE_coordinates);

    % Number of UEs
    numUEs = K; 

    % Number of clusters
    newClusterSize = ceil(numUEs / pilotLength);

    % Initial centroids
    centroids = datasample(data, newClusterSize, 'Replace', false);
    disp(centroids);
    
    % Initialize convergence flag
    converged = false;

    % Set convergence threshold
    convergenceThreshold = 0.001;

    % Initialize iteration counter
    iterations = 0;

    % Initialize cluster indices
    clusterIndices = zeros(numUEs, 1);

    % Set maximum number of iterations
    maxIterations = 1000;
    clustershrink = false;

%********************************************************************************
%% *Cluster formation using K-means and handling UEs in each cluster dynamically*
%********************************************************************************
    
    % Calculate distances from centroids
    distances = pdist2(data, centroids);
    
    % Assign points to clusters based on minimum distance
    [~, clusterIndices] = min(distances, [], 2);
    
    %% Continue iterating until convergence
    while ~converged
        % Increment iteration count for performance check
        iterations = iterations + 1;
        
        % Update centroids
        oldCentroids = centroids;
        for i = 1:newClusterSize
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

      %% Check for clusters exceeding 10 UEs and increasing cluster size OR 
      %% skip the process in case iterations goes beyond 1000 and form clusters

      if converged && clustershrink
         clustershrink = false; % Reset the clustershrink flag
      else
         if converged
           clusterSizes = zeros(1, newClusterSize);
           for i = 1:newClusterSize
               clusterSizes(i) = sum(clusterIndices == i);
           end
           largeClusters = find(clusterSizes > pilotLength);
           if ~isempty(largeClusters)
              % Increase new cluster size count by one
              newClusterSize = newClusterSize + 1; % Increment by a fixed value
              centroids = datasample(data, newClusterSize, 'Replace', false);
              % Update centroids
              oldCentroids = centroids;

              % Calculate distances from centroids
              distances = pdist2(data, centroids);
            
              % Assign points to clusters based on minimum distance
              [~, clusterIndices] = min(distances, [], 2);
    
              % Reset convergence flag to continue iterations
              converged = false;
           end
         end
      end
        
      %% If convergence is not met within maximum iterations, resetting flags to skip
      %% the process of forming clusters exceeding 10 UEs

      if ~converged && iterations >= maxIterations
         if newClusterSize > 1
           newClusterSize = newClusterSize - 1; % Shrink cluster size
           centroids = datasample(data, newClusterSize, 'Replace', false);

           % Update centroids and distances
           oldCentroids = centroids;

           % Calculate distances from centroids
           distances = pdist2(data, centroids);

           % Assign points to clusters based on minimum distance
           [~, clusterIndices] = min(distances, [], 2);

           % Reset flags
           converged = false;
           iterations = 0;
           clustershrink = true;
         end
      end

      % Place this initialization block before its usage
      finalClusters = cell(newClusterSize, 1);
      for i = 1:newClusterSize
          finalClusters{i} = find(clusterIndices == i);
      end
      %% Identify empty clusters
      emptyClusters = find(cellfun(@isempty, finalClusters));
        
      %% Remove empty clusters and adjust newClusterSize
      if ~isempty(emptyClusters)
          centroids(emptyClusters, :) = [];
          newClusterSize = newClusterSize - numel(emptyClusters);
          distances = pdist2(data, centroids);
          [~, clusterIndices] = min(distances, [], 2);
      end
    end
    
%**************************************************************************
%% **************Initialize cell array to store cluster indices************
%**************************************************************************
    clusterUEIndices = cell(newClusterSize, 1);

    % Store cluster indices
    for i = 1:newClusterSize
        clusterUEIndices{i} = find(clusterIndices == i);
    end
    
    % Display cluster indices
    for i = 1:numel(clusterUEIndices)
        fprintf('ClusterIndices %d:\n', i);
        disp(clusterUEIndices{i});
    end

    cluster_matrix = zeros(pilotLength, newClusterSize);
    % Display the final cluster indices
    for i = 1:newClusterSize
        clusterUEs = find(clusterIndices == i);
        %cluster_matrix(:, i) = find(clusterIndices == i);
        cluster_matrix(1:numel(clusterUEs), i) = clusterUEs;
        fprintf('Indices in Cluster %d: %s\n', i, num2str(clusterUEs'));
    end

%**************************************************************************
%% ************************Pilot Assignment********************************
%**************************************************************************
extractedPilotIndices = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C,cluster_matrix,1);
