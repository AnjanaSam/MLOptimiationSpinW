% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function option = InitializeParticleSwarm_FcutoffSR(ub, lb, nPara, numParticles, Pcutoff,Plevel,fac,SReg)

Iteration = 0;
% Pcutoff = 0.5;


%%...Alogithm Starts here.................................................
%.........................................................................

lbMatrix = repmat(lb,numParticles,1);
ubMatrix = repmat(ub,numParticles,1);

SelectedRegion = SReg;

ParamSR = SelectedRegion(:,1:end-1);
SelSM = SelectedRegion(:,end);

if size(SReg,1)>1
    
    
    
    sParam = ParamSR(SelSM==0,:);
    mParam = min(sParam, [],1);
    MParam = max(sParam, [], 1);
    dParam = 1.5*(MParam - mParam)/2;
    ParamMean = (MParam + mParam)/2;
    
    ubs = ParamMean + dParam;
    lbs = ParamMean - dParam;
    
    vmax = (ub-lb);
    
    ub(vmax==0) = ubs(vmax==0);
    lb(vmax==0) = lbs(vmax==0);
    
    ub = min([ub;ubs]);
    lb = max([lb;lbs]);
end


% Evaluated Function Values
vmax = (ub-lb);

Position = zeros(numParticles,nPara);
con = 1;
iT = 0;
mm = 1;

while con
    
    pos = repmat(lb,fac*numParticles,1) + rand(fac*numParticles,nPara).*repmat(vmax,fac*numParticles,1);
    
    [IdxSM, ~] = knnsearch(ParamSR(:,vmax~=0),pos(:,vmax~=0));
    selsm = SelSM(IdxSM);
    
    
    accptance = (selsm==0);
    
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
    
    mm = mm+1;
    
    
end


% Positions = repmat(lb,numParticles,1) + rand(numParticles,nPara).*repmat(vmax,numParticles,1);


option = struct('Iteration',Iteration,'ub',ub, 'lb',lb, 'nPara',nPara, 'numParticles',numParticles, ...
    'lbMatrix',lbMatrix,'ubMatrix', ubMatrix,...
    'Positions',Position,'AllSample',[],'AllFun',[],'Pcutoff',Pcutoff,'Plevel',Plevel,'fac',fac,'SelectedRegion',SReg);

