% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function [XX,YY,Z, InsI] = readSimDataSym(SFolder,I0, datafile, Range)

% datafile = '../Na2NiTeO6_5K_15meV_bkg.bin';
% datafile = '../Na2NiTeO6_5K_15meV_bkg2.bin';

Data = readBin(datafile);
lv = Data(1);
sA = reshape(Data(1+[1:2*lv]),2,lv);
psA = prod(sA);

is = 2*lv+1;
q = reshape(Data(is+[1:psA(1)]),sA(1,1),sA(2,1))';
is = is+psA(1);
de = reshape(Data(is+[1:psA(2)]),sA(1,2),sA(2,2))';
is = is+psA(2);
ZZ = reshape(Data(is+[1:psA(3)]),sA(1,3),sA(2,3));

Ind = I0;
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

InsI = Ins;

Z = ZZ;

[XX,YY] = meshgrid(K,Omega);

Mask = (YY<Range(1))|(YY>Range(2));
Z(Mask) = nan;
%         bZ = Z(0.9<YY&YY<1);
bZ = 0;
Z = Z - mean(bZ(~isnan(bZ)));

MaskG = ~isnan(Z);

nDataD = std(Z(MaskG));
Z = Z / nDataD;



AAlpha = 1;
InsI(Mask) = nan;
nData = std(InsI(MaskG&(~isnan(InsI))));
InsI = InsI/nData;















