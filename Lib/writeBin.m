% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function writeBin(paramName,para0V)

lpV = length(para0V);
sA = zeros(2,lpV);
for i1 = 1:lpV
    val = para0V{i1};
    
    vatt = size(val);
    sA(:,i1) = vatt(1:2);
end

if coder.target('MATLAB')
    fileID = fopen(paramName,'w');
    fwrite(fileID,lpV,'float64');
    fwrite(fileID,sA,'float64');
    for i1 = 1:lpV
        fwrite(fileID,para0V{i1},'float64');
    end
    fclose(fileID);
else
    fileID = coder.opaque('FILE *', 'NULL','HeaderFile','<stdio.h>');
    fileID = coder.ceval('fopen', [paramName, char(0)], ['wb+', char(0)]);
    val = lpV; tVal = val([1:end]); len = length(tVal); coder.ceval('fwrite', coder.ref(tVal), 8,len,fileID);
    val = sA; tVal0 = val([1:end]); len = length(tVal0); coder.ceval('fwrite', coder.ref(tVal0), 8,len,fileID);
    mul = prod(sA,1);
    for i1 = 1:lpV
        if mul(i1)~=0
            tVal2 = zeros(mul(i1),1);
            tVal2(:) = reshape(para0V{i1},mul(i1),1); 
%             tVal2 = para0V{i1}(1:end); 
            coder.ceval('fwrite', coder.ref(tVal2), 8,mul(i1),fileID);
        end
    end
    coder.ceval('fclose',fileID);
end