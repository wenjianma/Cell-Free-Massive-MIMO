%simple Pilot Assignment algorithm
%input:
%   input1 == UE_table
%   input2 == total_pilots
%output: the pilot_index (1d array)

%IMPORTANT: give me the input of the 1d arrays as (1, i), NOT (i, 1)
%IMPORTANT2: output is given in the format of (1, i) as well

function [return_value1] = simple_alg(input1, input2)
    Q = input1;

    total_UE = width(Q);

    total_pilots = input2;

    pilot_index = zeros(1, total_pilots);
    
    %assign all the orthogonal pilots
    %checking for pilot_index as well, in case we have less users than
    %pilots
    i = 1;
    while i <= total_pilots && i <= width(pilot_index)
        pilot_index(i) = i;
        i = i + 1;
    end

    %reassign the pilots based on additional interference introduced by
    %each new assignment (pick the lowest amount)

    % i <= total_UE, because we have already assigned pilots to
    %the first users and i starts from total_pilots+1

    %important: this assumes that the pilot_index array is initialized with
    %zeros! (to stop the while once it hits the first 0 entry)
    i = total_pilots + 1;
    while i <= total_UE
        j = 2;
        %set this to a very big value
        pilot_index(i) = pilot_index(1);
        minimum_interference = calculate_interference(pilot_index);
        while j <= total_pilots && pilot_index(j)~=0
            
            %this is why it is a good idea to have the object store all the
            %values needed for the calculation, we would just pass the
            %object here...

            %to do: find out what input is needed here:
            additional_interference = calculate_interference(pilot_index);
            %fprintf('%d %f', additional_interference);
            if(additional_interference < minimum_interference)
                minimum_interference = additional_interference;
                pilot_index(i) = pilot_index(j);
            end
            j = j + 1;
        end
        i = i + 1;
    end

    %return the pilot_index
    return_value1 = pilot_index;
end

function [return_value1] = calculate_interference (input1, input2)
    %do stuff to calculate the additional interference
    %would love to have the objects as input
    return_value1 = randsample(1000, 1);
end
