% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function data = readBin(fileName,silent)

if coder.target('MATLAB')
    tfileID = fopen(fileName,'r');
    file_length = get_file_length(tfileID);
    data = fread(tfileID,[file_length,1],'float64');
    fclose(tfileID);
else
    tfileID = coder.opaque('FILE *', 'NULL','HeaderFile','<stdio.h>');
    tfileID = coder.ceval('fopen', [fileName, char(0)], ['r', char(0)]);
    
    lSize = coder.opaque('long', 'NULL');
    lSize0 = coder.opaque('long', 'NULL');
    
    coder.ceval('fseek',tfileID, 0, 2);
    lSize = coder.ceval('ftell', tfileID);
    coder.ceval('rewind', tfileID);
    lSize0 = coder.opaque('long','lSize/8');
    lSize1 = int32(lSize0);
    
    if silent
        coder.ceval('printf', ['%s', char(10), char(0)],[fileName, char(0)]);
        coder.ceval('printf', ['%d %d', char(10), char(0)],lSize,lSize0);
    end

    data = double(zeros(lSize1,1));coder.ceval('fread', coder.ref(data), 8,lSize1, tfileID); 
    coder.ceval('fclose',tfileID);
end
