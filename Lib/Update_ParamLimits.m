% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function MappingObj = Update_ParamLimits(MappingObj,ID,scalefac)

if nargin < 3
   scalefac = 10; 
end

AllSample = MappingObj.AllSample;
AllFun0 = MappingObj.AllFun;

ub = MappingObj.ub;
lb = MappingObj.lb;

Pcutoff = MappingObj.Pcutoff;

minF0 = min(AllFun0, [], 1);
maxF0 = max(AllFun0, [], 1);
dF0 = repmat(maxF0 - minF0,size(AllFun0,1),1);
mF0 = repmat(minF0,size(AllFun0,1),1);

AllFun0 = (AllFun0 - mF0)./dF0;
AllFun = sum(AllFun0(:,ID),2);
%     AllFun = AllFun0(:,1);
minF = min(AllFun);
maxF = max(AllFun);
AllFun = (AllFun - minF)/(maxF-minF);

SelectedRegion = MappingObj.SelectedRegion;

ParamSR = SelectedRegion(:,1:end-1);
SelSM = SelectedRegion(:,end);

vmax = (ub-lb);
[IdxSM, ~] = knnsearch(ParamSR(:,vmax~=0),AllSample(:,vmax~=0));
    selsm = SelSM(IdxSM);

% funN =   (AllFun- minF);

sel = (AllFun < scalefac*Pcutoff*(maxF - minF))&(selsm==0);

if sum(sel)>1
    SP = AllSample(sel,:);
    
    MappingObj.ub = min([max(SP,[],1);ub]);
    MappingObj.lb = max([min(SP,[],1);lb]);
end
