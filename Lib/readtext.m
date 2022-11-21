% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function fstr = readtext(file_name)

if coder.target('MATLAB')
    fstr = fileread(file_name);
else
    lfile = int16(0);
    fout = char(zeros(1,100));
    fout(1:length(file_name)) = file_name;
    
    fileID = coder.opaque('FILE *', 'NULL','HeaderFile','<stdio.h>');
    fileID = coder.ceval('fopen', [fout, char(0)], ['rb', char(0)]);
    
    % get the File Length..........
    coder.ceval('fseek', fileID, int32(0), coder.opaque('int', 'SEEK_END'));
    lfile = int32(0);
    lfile = coder.ceval('ftell', fileID);
    % Reset current file position
    coder.ceval('fseek', fileID, int32(0), coder.opaque('int', 'SEEK_SET'));
    
    
    % Remaining is the number of bytes to read (from the file)
    remaining = lfile;
    
    % Initialize a buffer
    buffer = zeros(1,lfile+1,'uint8');
    
    % Index is the current position to read into the buffer
    index = int32(1);
    
    while remaining > 0
        % Buffer overflow?
        if remaining + index > size(buffer,2)
            %             fprintf('Attempt to read file which is bigger than internal buffer.\n');
            %             fprintf('Current buffer size is %d bytes and file size is %d bytes.\n', size(buffer,2), filelen);
            break
        end
        % Read as much as possible from the file into internal buffer
        nread = coder.opaque('size_t');
        nread = coder.ceval('fread', coder.ref(buffer(index)), int32(1), remaining, fileID);
        n = int32(0);
        n = coder.ceval('(int)',nread);
        if n == 0
            % Nothing more to read
            break;
        end
        % Did something went wrong when reading?
        if n < 0
            %             fprintf('Could not read from file: %d.\n', n);
            break;
        end
        % Update state variables
        remaining = remaining - n;
        index = index + n;
    end
    
    % Close file
    coder.ceval('fclose', fileID);
    
    fstr = char(buffer(1:index));
end
end