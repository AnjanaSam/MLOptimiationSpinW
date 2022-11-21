SFolder = './Test/';
datafile = './kniaso4_13mev_300hz_5K_Bkg.bin';
% Energy Range
Range = [0.5,11];
% Order: K, J, G, Gp, J2, J3, Jz, A0
Param = [0   -0.7316         0         0   -0.1524    1.5560    0.6457   -0.3489]; 
Ind = 1; NQ = 300; 
EvalFlag = 0; % 1 -  Perform the powder averaged calculation



if ~exist(SFolder)
        mkdir(SFolder);
end

%_________________________________________________
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
%___________________________________________________
    

% Modify "calculate_sw2_RealS02_KNiAsO4_Opt04" script as needed
calculate_sw2_RealS02_KNiAsO4_Opt04(Param,Ind,SFolder,q,de, NQ, EvalFlag);

[XX,YY,Z, InsI] = readSimDataSym(SFolder,Ind, datafile, Range);

figure
subplot(1,2,1)
surf(XX,YY,Z,'edgecolor','none')
view(2)
setFigStandards
caxis([0.,0.6]*1e1)
grid off
% ylim([0,0.9])
% xlim([0,1.6])
title('Data')
set(gca,'fontsize',12,'linewidth',1)

subplot(1,2,2)
% InsI = imgaussfilt(InsI, 1);
InsI(isnan(Z)) = nan;
surf(XX,YY,InsI,'edgecolor','none')
view(2)
setFigStandards
caxis([0.,0.6]*1e1)
grid off
% ylim([0,0.9])
% xlim([0,1.6])
title('Sim')
set(gca,'fontsize',12,'linewidth',1)




