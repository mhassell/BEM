% a script to dump the quadrature rules into a text file in C++ format with
% std::vectors and shit.

% last modified: June 24, 2017

file = 'gaussianQuad.txt';
cases = [1:2:33 39 63];
fileID = fopen(file,'w');

for j = cases
    rule = tableGauss(j);
    fprintf(fileID, 'qd = {');
    for k = 1:size(rule,1)-1
        try
            fprintf(fileID, '{ %d, %d},\n', rule(k,1), rule(k,2));
        catch
            keyboard
        end
    end
    fprintf(fileID,'{ %d, %d}};\n\n', rule(end,1), rule(end,2));
end

fclose(fileID);

%% 

file = 'LogGaussianQuad.txt';
cases = [1:40];
fileID = fopen(file,'w');

for j = cases
    rule = tableLogGauss(j);
    fprintf(fileID, 'qd = {');
    for k = 1:size(rule,1)-1
        try
            fprintf(fileID, '{ %d, %d},\n', rule(k,1), rule(k,2));
        catch
            keyboard
        end
    end
    fprintf(fileID,'{ %d, %d}};\n\n', rule(end,1), rule(end,2));
end

fclose(fileID);