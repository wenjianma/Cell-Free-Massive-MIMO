%% General Description 
% This function is about generating the predetermined_clustering and the UE-Cluster association
% The UE-Clustering association is based on the aggregation of the channel gain over noise
% Basically, we divide each row of AP into one cluster
%% Function Description
function sortedMaxIndices = predetermined_clustering(K,L,bklT)
% Generate the predetermined clusters based on the number of AP
cluster_size = sqrt(L);
cluster_number = sqrt(L);

% for each UE, assign one cluster based on the bklT
for ue = 1:K
    ue_cluster = 0; maxbkl = -100000;
    % Go through all clusters
    for cluster_no = 1:cluster_number
        bklsum = 0;
        for kk = 1:cluster_size
            bklsum = bklsum + bklT(ue,(cluster_no-1)*cluster_size + kk);
        end
        
        if bklsum > maxbkl
            ue_cluster = cluster_no;
            maxbkl = bklsum;
        end
    end
    
    % Assign UE-Cluster association
    for kk1 = 1:cluster_size
        sortedMaxIndices(ue,kk1) = (ue_cluster-1)*cluster_size + kk1;
    end
end
    sortedMaxIndices = sortedMaxIndices';
end