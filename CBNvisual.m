function CBNvisual(input, output, Lx, Ly, start, incr, endVal)
    % Function to convert the provided Java code to MATLAB
    if nargin < 7
        start = 0;
        incr = 1;
        endVal = 0;
    end
    
    % Loop over the range
    for m = start:incr:endVal
        % Modify input and output file names based on the current index (m)
        if contains(input, '#')
            input = strrep(input, '#', num2str(m));
        else
            input = strcat(input, num2str(m), '.txt');
        end
        
        if contains(output, '#')
            output = strrep(output, '#', num2str(m));
        else
            output = strcat(output, num2str(m), '.png');
        end
        
        if ~exist(input, 'file')
            disp(['File ' input ' not found!']);
            return;
        else
            % Read the data from the input file
            fid = fopen(input, 'r');
            line = fgetl(fid);
            parts = str2double(strsplit(line));
            if length(parts) ~= 6
                disp('Invalid data!');
                fclose(fid);
                return;
            end
            
            % Initialize the data array
            data = zeros(Ly, Lx, 6);
            
            % Variables for tracking min/max values
            min_C = 1.0E100;
            max_C = -1.0E100;
            min_B = 1.0E100;
            max_B = -1.0E100;
            min_N = 1.0E100;
            max_N = -1.0E100;
            min_u_C = 1.0E100;
            max_u_C = -1.0E100;
            min_u_BN = 1.0E100;
            max_u_BN = -1.0E100;
            
            % Read the rest of the data
            for j = 1:Ly
                for r = 1:Lx
                    line = fgetl(fid);
                    temp = str2double(strsplit(line));
                    C = temp(1);
                    B = temp(2);
                    N = temp(3);
                    u_C = temp(4);
                    u_BN = temp(5);
                    data(j, r, :) = [C, B, N, u_C, u_BN, temp(6)];
                    
                    % Update min/max values
                    min_C = min(min_C, C);
                    max_C = max(max_C, C);
                    min_B = min(min_B, B);
                    max_B = max(max_B, B);
                    min_N = min(min_N, N);
                    max_N = max(max_N, N);
                    min_u_C = min(min_u_C, u_C);
                    max_u_C = max(max_u_C, u_C);
                    min_u_BN = min(min_u_BN, u_BN);
                    max_u_BN = max(max_u_BN, u_BN);
                end
            end
            
            fclose(fid);
            
            % Create the image
            image = zeros(Ly, Lx, 3, 'uint8');
            
            for j = 1:Ly
                for i = 1:Lx
                    C = data(j, i, 1);
                    B = data(j, i, 2);
                    N = data(j, i, 3);
                    u_C = data(j, i, 4);
                    u_BN = data(j, i, 5);
                    
                    r = uint8(255 * map(u_BN, min_u_BN, max_u_BN) * map(B, min_B, max_B));
                    g = uint8(255 * map(u_BN, min_u_BN, max_u_BN) * map(N, min_N, max_N));
                    b = uint8(255 * map(u_C, max_u_C, min_u_C) * map(C, min_C, max_C));
                    
                    image(Ly - j + 1, i, :) = [r, g, b];  % Flip the y-axis to match Java behavior
                end
            end
            
            % Save the image as a PNG
            imwrite(image, output);
        end
    end
end

% Helper function to map values between 0 and 1
function result = map(x, minVal, maxVal)
    result = (x - minVal) / (maxVal - minVal);
end
% to run the code: 
%CBNvisual('(txt data file)', 'outputFile.png', 1024, 1024, 0, 1, 10);
%
