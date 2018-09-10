function compressed_binary_inputs = compress_inputs(binary_inputs)
% compresses binary inputs into a starting state (1 or 0) and a list of time
% points at which the states switch into an opposite state 
N = size(binary_inputs, 1);
compressed_binary_inputs = struct([]);
for j = 1:N
    compressed_binary_inputs(j).start = binary_inputs(j,1);
end
for j = 1:N
    compressed_binary_inputs(j).switches = [];
    for i = 2:size(binary_inputs, 2) 
        if abs(binary_inputs(j,i) - binary_inputs(j,i-1)) > 1e-6
            compressed_binary_inputs(j).switches = [compressed_binary_inputs(j).switches, i];
        end
    end
end