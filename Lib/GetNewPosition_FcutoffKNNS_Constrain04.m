% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function MappingObj = GetNewPosition_FcutoffKNNS_Constrain04(MappingObj,myConstrain,BoundaryFun,ID)

if nargin<2
    myConstrain = [];
end

if nargin<3
    BoundaryFun = [];
end


% Iteration = MappingObj.Iteration;
ub = MappingObj.ub;
lb = MappingObj.lb;

fac = MappingObj.fac;

vmax = (ub-lb);

% if Iteration == 0
%
% else

numParticles = MappingObj.numParticles;
AllSample = MappingObj.AllSample;
SelectedRegion = MappingObj.SelectedRegion;

ParamSR = SelectedRegion(:,1:end-1);
SelSM = SelectedRegion(:,end);

D0 = 0.0025*sqrt(sum(vmax.^2));

AllFun0 = MappingObj.AllFun;

if nargin<4
    ID = 1:size(AllFun0,2);
end

Pcutoff = MappingObj.Pcutoff;
nPara = MappingObj.nPara;
con = 1;
Position = zeros(numParticles,nPara);
iT = 0;
m = 1;

if ~isempty(AllSample)
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
    cutoff = Pcutoff*(maxF - minF);
else
    cutoff = 1;
end

mm = 0;

while con
    
    if ~isempty(BoundaryFun)
        pos = BoundaryFun(fac*numParticles, lb, ub);
    else
        pos = repmat(lb,fac*numParticles,1) + rand(fac*numParticles,nPara).*repmat(vmax,fac*numParticles,1);
    end
    
    if ~isempty(AllSample)
        [Idx, D] = knnsearch(AllSample(:,vmax~=0),pos(:,vmax~=0));
        funN = AllFun(Idx);
    end
    
    [IdxSM, ~] = knnsearch(ParamSR(:,vmax~=0),pos(:,vmax~=0));
    selsm = SelSM(IdxSM);
    
    
    %         funN = griddatan(AllSample(:,vmax~=0),AllFun,pos(:,vmax~=0),'nearest');
    %         [funN,~,~] = predict(gprMd0,pos(:,vmax~=0));
    
    %         sel = ((funN- minF)>cutoff)&(~isnan(funN))&(D<D0); % Correct
    if ~isempty(AllSample)
        sel = (isnan(funN))|(D>D0);
        %         accptance = sel & (exp(1-funN/cutoff) > rand(fac*numParticles,1)) & (selsm==0);
        % Sample more around the Boundary
        accptance = sel & (exp(1-abs(cutoff-funN)/0.25/cutoff) > rand(fac*numParticles,1)) & (selsm==0);
    else
        accptance = (selsm==0);
        
    end
    if ~isempty(myConstrain)
        selcon = myConstrain(pos(accptance,vmax~=0));
        accptance(accptance) = selcon;
    end
    
    nsel = sum(accptance);
    
    if nsel==0
        mm = mm +1;
    else
        mm = 0;
    end
    
    if mm == 100
        warning('Idel interation limit reached (by Anjana)\r')
        
        pos2 = pos;
        nsel = fac*numParticles;
    else
        pos2 = pos(accptance,:);
    end
    
    iT0 = iT;
    if (numParticles - iT)<=nsel
        iLim = (numParticles - iT);
        iT = numParticles;
    else
        iLim = nsel;
        iT = iT +nsel;
    end
    
    Position(iT0+[1:iLim],:) = pos2([1:iLim],:);
    
    if iT==numParticles
        con  =0;
    end
    fprintf('Itr = %d Filled %d/%d\r',m,iT,numParticles);
    m = m+1;
end

%     funN = griddatan(AllSample(:,vmax~=0),AllFun,Position(:,vmax~=0),'nearest');
%     sel = ((funN- minF)>cutoff)&(~isnan(funN));

%     if sum(sel)~=0
%         error('Anjana')
%     else

%     end

MappingObj.Positions = Position;

% figure;
% scatter3(AllSample(:,3),AllSample(:,4),AllSample(:,2),[],AllFun,'filled'); hold on
% II = Position;plot3(II(:,3),II(:,4),II(:,2),'kx')
% xlabel('J_3')
% ylabel('J_3''')
% zlabel('J_2')
% rotate3d on
% xlim([-1,1]*0.3)
% ylim([0,1]*0.3)
% zlim([-1,1]*0.6)
% 
% if ~isempty(AllSample)
%     fprintf('[%.3f, %.3f] %.3f, m = %d\r',maxF,minF,cutoff,m-1)
end


% end






















