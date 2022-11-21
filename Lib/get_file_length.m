% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function file_length = get_file_length(fid)
% extracts file length in bytes from a file opened by fopen
% fid is file handle returned from fopen

% store current seek
current_seek = ftell(fid);
% move to end
fseek(fid, 0, 1);
% read end position
file_length = ftell(fid);
% move to previous position
fseek(fid, current_seek, -1);

end