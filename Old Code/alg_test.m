
%create a UE array:

%first value in input is number of rows, in this case, 10
%second value is number of columns
%if second value is skipped, create 2d with equal sides

%objBArray = repmat(objB(), size(objAArray));
UE_array = repmat(UE_packet(), 1, 10);

i = 1;

while i <= 10
    j = 1;
    while j <= 2
        UE_array(i, j) = UE_packet(i+j, 10, 5);
        j = j + 1;
    end
    i = i + 1;
end

%pilot_array = zeros(1, 10);

pilot_array = simple_alg(UE_array, 5);
i = 1;
while i<width(pilot_array)
    fprintf('%d %f', pilot_array(1, i));
    i = i + 1;
end
