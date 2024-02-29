function movingSum = movingCumulSum(inputArray, windowSize)
% Check if the window size is valid
if windowSize <= 0 || windowSize > length(inputArray)
    error('Window size must be a positive integer smaller than or equal to the length of the input array.');
end

conv_factor = exp(-[0:1/24:1e4]/windowSize);
conv_factor = conv_factor / sum(conv_factor);
% Initialize the output array to store the moving cumulative sum
movingSum = conv(inputArray,conv_factor,'full');
movingSum = movingSum(1:length(inputArray));
end