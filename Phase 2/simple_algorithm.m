%simple Pilot Assignment algorithm
%input:
%   input1 == UE_table
%   input2 == total_pilots
%output: the pilot_index (1d array)

%IMPORTANT: give me the input of the 1d arrays as (1, i), NOT (i, 1)
%IMPORTANT2: output is given in the format of (1, i) as well

% function [return_value1] = simple_alg(input1, input2)
%     Q = input1;
%
%     total_UE = width(Q);
%
%     total_pilots = input2;
%
%     pilot_index = zeros(1, total_pilots);
%
%     %assign all the orthogonal pilots
%     %checking for pilot_index as well, in case we have less users than
%     %pilots
%     i = 1;
%     while i <= total_pilots && i <= width(pilot_index)
%         pilot_index(i) = i;
%         i = i + 1;
%     end
%
%     %reassign the pilots based on additional interference introduced by
%     %each new assignment (pick the lowest amount)
%
%     % i <= total_UE, because we have already assigned pilots to
%     %the first users and i starts from total_pilots+1
%
%     %important: this assumes that the pilot_index array is initialized with
%     %zeros! (to stop the while once it hits the first 0 entry)
%     i = total_pilots + 1;
%     while i <= total_UE
%         j = 2;
%         %set this to a very big value
%         pilot_index(i) = pilot_index(1);
%         minimum_interference = calculate_interference(pilot_index);
%         while j <= total_pilots && pilot_index(j)~=0
%
%             %this is why it is a good idea to have the object store all the
%             %values needed for the calculation, we would just pass the
%             %object here...
%
%             %to do: find out what input is needed here:
%             additional_interference = calculate_interference(pilot_index);
%             %fprintf('%d %f', additional_interference);
%             if(additional_interference < minimum_interference)
%                 minimum_interference = additional_interference;
%                 pilot_index(i) = pilot_index(j);
%             end
%             j = j + 1;
%         end
%         i = i + 1;
%     end
%
%     %return the pilot_index
%     return_value1 = pilot_index;
% end
%
% function [return_value1] = calculate_interference (input1, input2)
%     %do stuff to calculate the additional interference
%     %would love to have the objects as input
%     return_value1 = randsample(1000, 1);
% end



%simple Pilot Assignment algorithm
% This is the algorithm in Ref no.7
%% Method 1:
% Input:
% R -> spatial correlation matrix N*N*K*L
% pilotlength -> the number of the orthogonal pilots
% K
% L
% Output:
% pilotindex -> K*1 matrix about the pilot assignment for each UE
function pilotIndex = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C)
%initialization:
pilotIndex = zeros(K,1);
%if the number of UE is smaller or equal to the number of pilot
if  K <= pilotLength
    for i = 1:K
        pilotIndex(i) = i;
    end
    % if the number of UE is larger than the number of pilot
else
    % first assign orthogonal pilots to first pilotlength UEs
    for i = 1: pilotLength
        pilotIndex(i,1) = i;
    end
    % then assign other pilots
    for i = pilotLength+1 : K
        minimum_interference = 10000000; id = 0;

        % shared_minimum_interference is to be used in the worst case
        % scenario, where all pilots are already shared and we cannot
        % guarantee separation of service
        % shared_pilot stores the pilot for that scenario
        % guard_3 keeps track of whether we have found a pilot that can
        % guarantee separation of clusters
        shared_minimum_interference = 10000000;
        shared_pilot = 0;
        non_shared_pilot = 0;

        % guard_1 keeps track of whether we have assigned a pilot that
        % is not being used by any other UE (maybe will remove this in
        % the future, not very useful if we assign the first orthogonal
        % pilots sequentially)
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
            y = 1;
            while y < length(pilotIndex) + 1
                % disp(y);
                % disp(pilotIndex(y, 1));
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
                % if the UEs that will share the pilot are
                % served by at least one common AP, then
                % guard_2 = 1, which means DO NOT assign the
                % pilot

                % this if assumes a sorted table (vertically)
                % if there is overlap:
                if isempty(intersect(sortedMaxIndices(:, using_pilot(y)), sortedMaxIndices(:, i))) == 0 
                    %if sortedMaxIndices(1,using_pilot(y)) == sortedMaxIndices(1,i) | sortedMaxIndices(2, using_pilot(y)) == sortedMaxIndices(2,i) | sortedMaxIndices(1, using_pilot(y)) == sortedMaxIndices(2,i)
                    % here we need to reset the
                    % minimum_interference and make sure to get
                    % a non AP-sharing pilot; worst case
                    % scenario, we keep the shared pilot that
                    % introduces the least amount of
                    % interference and use that
                    
                    % if the overlap AP is not AP 0, we do put the overlap flag
                    overlapAPs = intersect(sortedMaxIndices(:, using_pilot(y)), sortedMaxIndices(:, i));
                    if sum(overlapAPs) > 0 
                    % pilot_counter counts how many pilots have overlap
                        pilot_overlap_flag = 1;
                    end
                    %% commented this section to fix the logical mistake discussed in the 20/11/2023 meeting

                    % if interference < shared_minimum_interference
                    %     shared_minimum_interference = interference;
                    %     shared_pilot = p;
                    % end

                    %% commented this section to fix the logical mistake discussed in the 20/11/2023 meeting
                    % else
                    %     % else, if there is no overlap:
                    %     if interference < minimum_interference
                    %         minimum_interference = interference;
                    %         non_shared_pilot = p;
                    %     end
                end

            end
            %% fixed the issue from 20/11/2023
            % first we go through all the UEs (that are sharing the same pilot)
            % to figure out if there is overlap,
            % and then proceed like this:
            if pilot_overlap_flag == 1
                pilot_counter = pilot_counter + 1;
                % for y = 1:using_pilot
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
                % end
            else
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
        %if pilot_counter >= length(using_pilot)
        if pilot_counter >= pilotLength
            id = shared_pilot;
        else
            id = non_shared_pilot;
        end
        pilotIndex(i,1) = id;
    end
end
end