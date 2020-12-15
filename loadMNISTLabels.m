function labels = loadMNISTLabels(filename)
%loadMNISTLabels returns a [number of MNIST images]x1 matrix containing
%the labels for the MNIST images
fp = fopen(filename, 'rb');
fseek(fp,8,'bof');
labels = fread(fp, inf, 'unsigned char');
fclose(fp);
end
%assert(fp ~= -1, ['Could not open ', filename, '']);
%magic = fread(fp, 1, 'int32', 0, 'ieee-be');
%assert(magic == 2049, ['Bad magic number in ', filename, '']);
%numLabels = fread(fp, 1, 'int32', 0, 'ieee-be');
%assert(size(labels,1) == numLabels, 'Mismatch in label count');