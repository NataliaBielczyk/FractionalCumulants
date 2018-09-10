function [binary_inputs1] = decompress_inputs(binary_inputs, N, T, Fs, high_input)
% decompresses binary inputs from compressed version into a~complete time series of 0s and 1s
binary_inputs1 = zeros(N, T*Fs);
for j = 1:N
    binary_inputs1(j,:) = repmat(binary_inputs(j).start, 1, T*Fs);
    actual = binary_inputs(j).start;
    sequence = binary_inputs(j).switches;
    for k = 1:length(sequence)
        binary_inputs1(j, sequence(k):end) = high_input - actual;
        actual = high_input - actual;
    end
end