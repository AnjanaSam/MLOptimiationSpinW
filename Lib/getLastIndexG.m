% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function Temp = getLastIndexG(SFolder,fileName)

file_name = [SFolder,'/',fileName];
fstr = readtext(file_name);

if isempty(fstr)
   Temp = 0;
else
    d = dlmread(file_name);
    Temp = d(end,1);
end


