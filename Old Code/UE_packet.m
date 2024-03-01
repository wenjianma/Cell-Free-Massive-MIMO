classdef UE_packet
    %UE_PACKET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pilot
        id
        coordinates_x
        coordinates_y
        separated_flag
        serving_APs
        UE_packet_k_dash
    end
    
    methods
        %this is the constructor 
        function obj = UE_packet(input1, input2, input3)
            if nargin>2
                %pilot = null;
                obj.separated_flag = false;

                %16 maximum serving_APs per UE
                obj.serving_APs = zeros(0, 16);

                obj.id = input1;
                obj.coordinates_x = input2;
                obj.coordinates_y = input3;
            
                %the rest of the inputs need to be the elements needed to
                %calculcate the SINR
            else
                %set everything to 0/false default constructor;
                obj.separated_flag = false;
                obj.serving_APs = zeros(16);
                obj.id = 0;
                obj.coordinates_x = 0;
                obj.coordinates_y = 0;
            end
        end
        

        %function for calculating the SINR of k, k'

        function return_value1 = calculate_SINR(input1, input2)
            %I honestly don't know why we need input2 for this to run

            %need to calculate SINR between two UE's (SINR<k, k'>)
            %we will treat .obj (the current object) as UE_packet_k
            %the input will be UE_packet_k_dash (k')
            obj.UE_packet_k_dash = input1;

            final_SINR_calculation = 1 + 2 + 3;

            return_value1 = final_SINR_calculation;
            
        end
    end
end

