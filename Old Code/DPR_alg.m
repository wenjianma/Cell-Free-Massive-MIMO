%BEFORE CALLING THIS ALGORITHM:
% 1. make sure that the total number of UE's are more than the available
% pilot sequences.
% 2. the created UE pairs must NOT be served by the same AP

function [return_value1] = DPR_alg(input1, input2, input3)
    
    %Input1 = Q table (2d array, containing all possible UE pairs
    % already sorted in descending order in terms of separation distance)
    %Keep in mind: each user in the UE pairs MUST be served by a different AP

    %maybe we can create a packet (class) called UE_packet and feed it into
    %this array. This packet should contain all of the information
    %necessary to calculate the SINR

    Q = input1;

    %Input2 = total number of UE's

    total_UE = input2;

    %initialize the Q table for testing purposes%
    %2 rows, x columns
    %x = 5;
    %Q = zeros(2, x);
    
    %Q2 table (2d array [total_UE/2], containing all the calculations of SINR<k, k'> &&
    %SINR<k', k>). This will have a maximum length equal to the total
    %number of UE's divided by 2, since we are comparing pairs

    %initialize the Q2 table%
    Q2 = zeros(2, total_UE/2);
    %Q2 = UE_packet.empty(0, total_UE/2);
    %Q2 = repmat(UE_packet(), 0, total_UE/2);

    %Q3 table (2d array [2 BY total_UE/2], containing all of the already selected pairs of users. This
    %will be used to filter out the already selected users and remove them
    %from table Q (to make searching faster and avoid confusion))

    %Maximum length should be equal to total_UE/2

    %Q3 = zeros(2, total_UE/2);
    %Q3 = UE_packet.empty(2, total_UE/2);
    Q3 = repmat(UE_packet(), 2, total_UE/2);

    %NOTE: Q2 and Q3 should have a 1-1 mapping of SINR results and UE pairs
    
    %Input3 = pilot_array (1d array [number_of_pilot_sequences] 
    %[otherwise called "vector"], where all
    %the pilot sequences are stored)

    pilot_array = input3;
    
    %return_value1 = assigned_pilots (1d array [total_UE] where we
    %store the UE_packet instances). The pilot can then be individually
    %retrieved using the get_pilot() method.
    %initialize it here, finalize it at the end of the function:

    %assigned_pilots = zeros(total_UE);
    %assigned_pilots = UE_packet.empty(0, total_UE);
    assigned_pilots = repmat(UE_packet(), 0, total_UE);

    %i = 1;
    %i is not needed, as the table has only 2 rows, so we can just go 1 and
    %then 1+1 for the indexing
    j = 1;
    %z is indexing the pilot_array, increasing by 1 every time we meet the
    %SINR requirement
    %z is also indexing the Q2 and Q3 arrays
    z = 1;

    %set the gamm_threshold to 5dB for testing, later we HAVE to calculate
    %this, either by using the IGS method, or my proposed binary search
    gamma_threshold = 5;

    %height(Q) returns the rows of the table - i
    %width(Q) returns the columns of the table - j

    %while i < height(Q)
        while j <= width(Q) & Q(j).coordinates_x ~= 0
            % SINR_1 = Q(1, j).calculate_SINR(Q(2, j));
            % SINR_2 = Q(2, j).calculate_SINR(Q(1, j));
            SINR_1 = Q(1, j).calculate_SINR(4);
            SINR_2 = Q(2, j).calculate_SINR(4);

            if SINR_1 > gamma_threshold && SINR_2 > gamma_threshold
                Q(1, j).pilot = pilot_array(z);
                Q(2, j).pilot = pilot_array(z);

                %store the result of the SINR calculcation in the Q2 table
                Q2(1, z) = SINR_1;
                Q2(2, z) = SINR_2;

                %store the pair in the Q3 table of selected pairs
                Q3(1, z) = Q(1, j);
                Q3(2, z) = Q(2, j);

                %remove the all instances of the pair from Q table
                %IMPORTANT: need to complete this part of the code to
                %remove all instances of these pairs from the Q table, pending...

                %maybe this will work:
                i1 = 1;
                while i1 < width(Q)
                    if Q(1, i1).id == Q(1, j).id || Q(1, i1).id == Q(2, j).id
                        Q(1, i1) = NaN;
                    else if  Q(2, i1).id == Q(1, j).id || Q(2, i1).id == Q(2, j).id
                        Q(2, i1) = NaN;
                    end
                    end
                    i1 = i1 + 1;
                end

                z = z + 1;
            end

            j = j + 1;
        end

    %end
    

    %if at least two users are left without a pilot, because they did not meet the
    %SINR threshold AND there are at least 2 unused pilots:
    %allocate the remaining pilots to the users
    %add them to a Q4 1d table
    %merge the Q3 and Q4 tables later
    %to be completed, pending...


    %if all pairs meet the SINR threshold (z == total_UE/2) AND  there is
    %at least 1 unused pilot
    %find the pairs in Q2 with the
    %lowest sum of SINR's and allocate a unique pilot to each of the users

    %need to index the pilot_array separately, because we use z for the
    %control statements (1st while)
    %basically, z reflects the amount of pilot sequences that we have
    %already allocated
    pilot_array_index = z;
    if z == total_UE/2 && width(pilot_array)-z >= 1
        counter = 1;
        %while we have more than one pilot left
        while counter <= width(pilot_array)-z
           j = 1;
           index = 1;
           min = Q2(1, j) + Q2(2,j);
           %find the pair with the lowest sum of SINR's and set the
           %separated_flag to true
           while j < width(Q2) && Q2.Value(j) ~= 0
               %check for min AND separated_flags to be false (~ == !)
                if ((Q2(1, j) + Q2(2,j)) < min) & ~(Q3(1, j).get_separated_flag()) & ~(Q3(2, j).get_separated_flag())
                    min = Q2(1, j) + Q2(2,j);
                    index = j;
                end
                j = j + 1;
           end
           %we can do this because Q2 and Q3 maintain a 1-1 index
           Q3(1, index).set_pilot(pilot_array.Value(pilot_array_index));
           pilot_array_index = pilot_array_index + 1;
           %set separated_flag to true for both UE's in pair:
           Q3(1, index).set_separated_flag(true);
           Q3(2, index).set_separated_flag(true);
           counter = counter + 1;
        end
    end

    %Now, turn the Q3 table into a 1d array and return it as
    %"assigned_pilots"
    
    i = 1;
    j = 1;

    while j <= width(Q3) & Q3(j).coordinates_x ~= 0
        assigned_pilots(i) = Q3(1, j);
        assigned_pilots(i+1) = Q3(2, j);

        i= i + 2;
        j= j + 1;
    end

    %TO BE DONE:
    %merge the assigned_pilots array with the Q4 table


    %return values finalized:
    return_value1 = assigned_pilots;
end
   
