% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................

% Log File Format..............
% Index Hamiltonian_Parameters Alpha Chi

function out = getChi2M_RealS_Sym_KNiAsO4R04(Param,SFolder,saveFlag, Range, NQ, datafile, smf, OutFile)


if nargin < 3
    saveFlag = 1;
    
end


LastI = getLastIndexG(SFolder,OutFile);

lP = size(Param,1);

NewInd = LastI+[1:lP]';
cOut = cell(lP,1);
IOut = cell(lP,1);

% p = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(p)
%     poolsize = feature('numcores');
%     p = parpool(poolsize);
% else
%     poolsize = p.NumWorkers;
% end
% fprintf('Num Cores = %d\n',poolsize);

% datafile = '../Na2NiTeO6_5K_15meV_bkg.bin';
Data = readBin(datafile);
lv = Data(1);
sA = reshape(Data(1+[1:2*lv]),2,lv);
psA = prod(sA);

is = 2*lv+1;
q = reshape(Data(is+[1:psA(1)]),sA(1,1),sA(2,1))';
is = is+psA(1);
de = reshape(Data(is+[1:psA(2)]),sA(1,2),sA(2,2));
is = is+psA(2);
ZZ = reshape(Data(is+[1:psA(3)]),sA(1,3),sA(2,3));



parfor i1 = 1:lP
%     for i1 = 1:lP
    
    
    [XX,YY] = meshgrid(q,de);
    
    param = Param(i1,:);
    Ind = NewInd(i1);
    
    tic
    
%     intensity_pwd = calculate_sw2_RealS02_KNiAsO4_Opt(param,Ind,SFolder, q, de, NQ);
%     intensity_pwd = calculate_sw2_RealS02_KNiAsO4_Opt02(param,Ind,SFolder, q, de, NQ);
    [intensity_pwd, magStrF] = calculate_sw2_RealS02_KNiAsO4_Opt04(param,Ind,SFolder, q, de, NQ);
    
    t = toc;
    fprintf('Porcess %d: Time = %.4f s \r\n',i1,t)
    
    NC = isnan(intensity_pwd);
    if sum(NC(:))==0
        Data = readBin([SFolder,'/Model_',num2str(Ind),'.bin']);
        
        lv = Data(1);
        
        sA = reshape(Data(1+[1:2*lv]),2,lv);
        psA = prod(sA);
        
        is = 2*lv+1;
        K = reshape(Data(is+[1:psA(1)]),sA(1,1),sA(2,1));
        is = is+psA(1);
        Omega = reshape(Data(is+[1:psA(2)]),sA(1,2),sA(2,2));
        is = is+psA(2);
        Ins = reshape(Data(is+[1:psA(3)]),sA(1,3),sA(2,3));
        
        %         [KK,OO] = meshgrid(K,Omega);
        
        %         InsI = interp2(KK,OO,Ins,XX,YY);
        InsI = Ins;
        if smf ~=0
            InsI = imgaussfilt(InsI, smf);
        end
        
        Z = ZZ;
        
        Mask = (YY<Range(1))|(YY>Range(2));
        Z(Mask) = nan;
        %         bZ = Z(0.9<YY&YY<1);
        bZ = 0;
        Z = Z - mean(bZ(~isnan(bZ)));
        
        MaskG = ~isnan(Z);
        
        nDataD = std(Z(MaskG));
        Z = Z / nDataD;
        
        
        
        AAlpha = 1;
        nData = std(InsI(MaskG&(~isnan(InsI))));
        InsI = InsI/nData;
        
        tt = (InsI-Z).^2;
        chi = mean(tt(:),'omitnan');
        
        %         bnds = [0.3,3,1];
        %
        %         [chi, AAlpha, ~, ~] = getMinChiPowder_go(XX,YY,InsI,Z,MaskG,bnds);
        %
        %
        %
        %         AA = roundn(AAlpha,-6);
        %         if (AA==bnds(1))|(AA==bnds(2))
        %             disp('Boundary Hit.....');
        %             eP = max(InsI,[],2);
        %             [~,locsS,~,~] = findpeaks(eP,1:length(Y),'Annotate','extents','MinPeakProminence',0.3 ...
        %                 ,'WidthReference','halfheight');
        %
        %             eD = max(Z,[],2);
        %             [~,locsD,~,~] = findpeaks(eD,1:length(Y),'Annotate','extents','MinPeakProminence',0.3 ...
        %                 ,'WidthReference','halfheight');
        %
        %
        %             if ((~isempty(locsS))&(~isempty(locsD)))
        %                 lb = Y(locsD(1))/Y(locsS(end));
        %                 ub = Y(locsD(end))/Y(locsS(1));
        %                 mb = Y(locsD(1))/Y(locsS(1));
        %
        %                 if (lb >= ub)|(lb<0.2)|(ub>5)
        %                     AAlpha = 1;
        %
        %                     nData = std(InsI(MaskG&(~isnan(InsI))));
        %                     InsI = InsI/nData;
        %
        %                     tt = (InsI-Z).^2;
        %                     chi = mean(tt(:),'omitnan');
        %
        %                 else
        %                     bnds = [lb,ub,mb];
        %                     [chi, AAlpha, ~, ~] = getMinChiPowder_go(XX,YY,InsI,Z,MaskG,bnds);
        %                 end
        %
        %             else
        %                 AAlpha = 1;
        %                 nData = std(InsI(MaskG&(~isnan(InsI))));
        %                 InsI = InsI/nData;
        %
        %                 tt = (InsI-Z).^2;
        %                 chi = mean(tt(:),'omitnan');
        %             end
        %         end
        
    else
        chi = 10;
        AAlpha = 1;
    end
    
    chiOri = chi;
    if ~magStrF
        chi = 10;
    
    end
    
    
    
    cOut{i1} = [chi, chiOri, AAlpha];
    %     sel = (0.24<XX&XX<0.32)&(0.61<YY&YY<0.67);
    %     IOut{i1} = sum(sum(Ins(sel)));
    
end



Out = cell2mat(cOut);
% Iou = cell2mat(IOut);
if saveFlag
    dlmwrite([SFolder,'/',OutFile],[NewInd,Param,Out],'-append')
end

out = Out(:,1);












































