
%% Method 1:
% Input:
% R -> spatial correlation matrix N*N*K*L
% pilotlength -> the number of the orthogonal pilots
% K
% L
% Output:
% pilotindex -> K*1 matrix about the pilot assignment for each UE
function pilotIndex = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C, cluster_matrix, alg_flag, extractedPilotIndices)

% alg_flag decides whether we are using k-means or not (k-means used,
% alg_flag == 1)

%initialization:
pilotIndex = zeros(K,1);
if alg_flag == 1
    pilotIndex = extractedPilotIndices;
end

%if the number of UE is smaller or equal to the number of pilot
if  K <= pilotLength
    for i = 1:K
        pilotIndex(i) = i;
    end
    % if the number of UE is larger than the number of pilot
else
    % % first assign orthogonal pilots to first pilotlength UEs
    % for i = 1: pilotLength
    %     pilotIndex(i,1) = i;
    % end
    % % then assign other pilots
    for i = 1:K
        minimum_interference = 10000000; id = 0;

        % shared_minimum_interference is to be used in the worst case
        % scenario, where all pilots are already shared and we cannot
        % guarantee separation of service
        % shared_pilot stores the pilot for that scenario
        shared_minimum_interference = 10000000;
        shared_pilot = 0;
        non_shared_pilot = 0;

        % pilot_counter counts how many pilots have AP overlap for a
        % specific User i
        pilot_counter = 0;

        for p = 1:pilotLength
            % calculate the interference defined in ref no.7
            % which is the contamination caused at the serving APs if
            % the p-th pilot is used by UE i
            interference = 0;

            %%%%
            % added for clustering
            % if the new minimal interference pilot is not used by
            % the UE's that share at least one common AP -> assign
            % this as the pilot
            % in the "using_pilot" array we store the UE's that are
            % sharing the same pilot, as the one we have just found
            % above
            using_pilot = zeros (K, 1);
            counter = 0;

            %for y = 1:pilotIndex
            % this for loop was broken for no reason, fuck Matlab...

            % This loop is used to find which UEs are using pilot p
            y = 1;
            while y < length(pilotIndex) + 1
                if p == pilotIndex(y, 1)
                    counter = counter + 1;
                    using_pilot(counter, 1) = y;
                end
                y = y + 1;
            end
            %%%%
            % extract serving APs from UEs (UEs are indeces to find their
            % serving APs from the sortedMaxIndices vector)
            % servingAPs dimensions [C, counter]
            servingAPs = sortedMaxIndices(:, using_pilot(1:counter));

            % Here we calculate the additional interference metric
            for kk = 1: counter
                for l = 1:C % with cluster formation, L will be changed (here L means the APs serving UE i)
                    if  ~(servingAPs(l,kk) == 0)
                        interference = interference + trace(R(:,:,using_pilot(kk),servingAPs(l,kk))*R(:,:,i,servingAPs(l,kk)));
                    end
                end
            end
            % update info


            % now that we have the array of the UEs that share this
            % pilot, we can find if these UEs are in the same (or
            % overlapping) clusters:
            pilot_overlap_flag = 0;
            for y = 1:using_pilot
                if using_pilot(y, 1) == 0
                    break
                end

                % if there is overlap:
                if isempty(intersect(sortedMaxIndices(:, using_pilot(y)), sortedMaxIndices(:, i))) == 0
                    %if sortedMaxIndices(1,using_pilot(y)) == sortedMaxIndices(1,i) | sortedMaxIndices(2, using_pilot(y)) == sortedMaxIndices(2,i) | sortedMaxIndices(1, using_pilot(y)) == sortedMaxIndices(2,i)

                    % if the overlap AP is not AP 0, we do put the overlap flag
                    % intersect returns all of the common values between
                    % the 2 arrays without repetition
                    overlapAPs = intersect(sortedMaxIndices(:, using_pilot(y)), sortedMaxIndices(:, i));
                    if sum(overlapAPs) > 0
                        pilot_overlap_flag = 1;
                    end
                end
            end
            %% fixed the issue from 20/11/2023
            % first we go through all the UEs (that are sharing the same pilot)
            % to figure out if there is overlap,
            % and then proceed like this:

            % if there is overlap for this particular pilot (after checking
            % which APs are serving all of the UEs that are also using this
            % pilot), check the shared_minimum_interference metric
            if pilot_overlap_flag == 1
                pilot_counter = pilot_counter + 1;
                if interference < shared_minimum_interference
                    shared_minimum_interference = interference;
                    shared_pilot = p;
                else
                    % if the metric is the same for both pilots, choose one
                    % randomly
                    if interference == shared_minimum_interference
                        zz = [p, shared_pilot];
                        pos = randi(length(zz));
                        shared_pilot = zz(pos);
                    end
                end

                % else, if there is no overlap for this pilot, check the
                % minimum_interference metric
            elseif pilot_overlap_flag == 0
                if interference < minimum_interference
                    minimum_interference = interference;
                    non_shared_pilot = p;
                else
                    % if the metric is the same for both pilots, choose one
                    % randomly
                    if interference == minimum_interference
                        zz = [p, non_shared_pilot];
                        pos = randi(length(zz));
                        non_shared_pilot = zz(pos);
                    end
                end
            end
        end
        % added for clustering
        % if there is still no pilot assigned, assign the one with the
        % least amount of interference (even if there is overlap)

        % for k-means | i == current UE
        % look for user i inside the cluster_matrix (cluster_matrix is
        % given as input from the k-means algorithm)
        % -> figure out in which cluster user i is placed
        [~,cluster_now] = find(cluster_matrix == i);
        swapdesh = 0;

        % this if is equivalent to the logical statement "if all of the
        % pilots have AP overlap, ..."
        if pilot_counter >= pilotLength
            % only used in k-means (alg_flag == 1)
            if alg_flag == 1

                % find which UEs are using the shared_pilot
                counter = 0;
                y = 1;
                using_pilot = zeros (K, 1);
                while y < length(pilotIndex) + 1
                    if shared_pilot == pilotIndex(y, 1)
                        counter = counter + 1;
                        using_pilot(counter, 1) = y;
                    end
                    y = y + 1;
                end

                for u = 1:using_pilot
                    % for all the users using the shared_pilot, find the one within
                    % the same k-means cluster as user i
                    % -> swap the pilots between the two users
                    [~, cluster_swap] = find(cluster_matrix == using_pilot(u));
                    if cluster_now == cluster_swap
                        temp = pilotIndex(using_pilot(u));
                        pilotIndex(using_pilot(u)) = pilotIndex(i);
                        pilotIndex(i) = temp;
                        swapdesh = 1;
                    end
                end
            end
            id = shared_pilot;

            % else, if at least one pilot does not have overlap
        else
            % for k-means | i == current UE
            % look for user i inside the cluster_matrix (cluster_matrix is
            % given as input from the k-means algorithm)
            % -> figure out in which cluster user i is placed
            % only used in k-means (alg_flag == 1)
            if alg_flag == 1
                % find which UEs are using the non_shared_pilot
                counter = 0;
                y = 1;
                using_pilot = zeros (K, 1);
                while y < length(pilotIndex) + 1
                    if shared_pilot == pilotIndex(y, 1)
                        counter = counter + 1;
                        using_pilot(counter, 1) = y;
                    end
                    y = y + 1;
                end
                for u = 1:using_pilot
                    % for all the users using the non_shared_pilot, find the one within
                    % the same k-means cluster as user i
                    % -> swap the pilots between the two users
                    [~, cluster_swap] = find(cluster_matrix == using_pilot(u));
                    if cluster_now == cluster_swap
                        temp = pilotIndex(using_pilot(u));
                        pilotIndex(using_pilot(u)) = pilotIndex(i);
                        pilotIndex(i) = temp;
                        swapdesh = 1;
                    end
                end
            end
            id = non_shared_pilot;

        end
        if swapdesh == 0
            pilotIndex(i,1) = id;
        end
    end
end
end