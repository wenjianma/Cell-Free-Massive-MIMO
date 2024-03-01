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
function pilotIndex = simple_algorithm(R,pilotLength,K,L,sortedMinIndices)
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
            guard = 0;
            for p = 1:pilotLength
                % calculate the interference defined in ref no.7
                % which is the contamination caused at the serving APs if
                % the p-th pilot is used by UE i
                interference = 0;
                for l = 1:L % with cluster formation, L will be changed (here L means the APs serving UE i)
                    ue_used = find(pilotIndex == p); % find the UEs which use the p-th pilot already
                    for kk = 1 : size(ue_used)
                        interference = interference + trace(R(:,:,ue_used(kk),l)*R(:,:,i,l));
                    end
                end
                % update info
                if interference < minimum_interference
                    minimum_interference = interference;
                    temp_var = p;
                    % added for clustering
                    % if the new minimal interference pilot is not used by
                    % the UE's that share at least one common AP -> assign
                    % this as the pilot
                    % in the "using_pilot" array we store the UE's that are
                    % sharing the same pilot, as the one we have just found
                    % above
                    using_pilot = zeros (K, 1);
                    counter = 0;
                    for y = 1:pilotIndex
                        if p == pilotIndex(y, 1)
                            counter = counter + 1;
                            using_pilot(counter, 1) = y;
                        end
                    end
                    % now that we have the array of the UEs that share this
                    % pilot, we can find if these UEs are in the same (or
                    % overlapping) clusters:
                    % guard_2 monitors if there is overlap in the clusters
                    % (at least one common AP)

                    guard_2 = 0;

                    % if counter > 0 ->  if the pilot is already used by at
                    % least one UE
                    if counter > 0
                        for y = 1:using_pilot
                            if using_pilot(y, 1) == 0
                                break
                            end
                            % if the UEs that will share the pilot are
                            % served by at least one common AP, then 
                            % guard_2 = 1, which means DO NOT assign the
                            % pilot

                            % this if assumes a sorted table (vertically)
                            if sortedMinIndices(1,using_pilot(y)) == sortedMinIndices(1,i) | sortedMinIndices(2, using_pilot(y)) == sortedMinIndices(2,i)
                                guard_2 = 1;
                            end
                        end
                    else
                        id = p;
                        guard = 1;
                    end
                    % if the clusters are not overlapping, assign the pilot
                    if guard_2 == 0
                        id = p;
                        guard = 1;
                    end
                end
            end
            % added for clustering
            % if there is still no pilot assigned, assign the one with the
            % least amount of interference (even if there is overlap)
            if guard == 0
                id = temp_var;
            end
            pilotIndex(i,1) = id;
        end
    end
end
