% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................

datafile = './kniaso4_13mev_300hz_5K_Bkg.bin';
SFolder = './DirectRunSym_kniaso4_1bK';
outSwam = 'DirectRun_RealS_Sym_kniaso4_1bK.mat';
Range = [0.5,11];

% Order: K, J, G, Gp, J2, J3, Jz, A0
PStr = {'K_1', 'J_1', '\Gamma_1', '\Gamma_1''', 'J_2', 'J_3', 'J_z', 'D'}


% This program can combine all the project folders and visulaize the
% results. Just add a new row for each project.
% Row: 'Project_folder', 1;
folders = {
    './DirectRunSym_kniaso4_1bK',1; % K-G
    };

cutoff = 1; % 0.04

%........................................................................
%.........................................................................

OutName= 'log_RealSMag.dat'; % cOut{i1} = [chi, chiOri, AAlpha];

D = [];
E = [];
SQW = {};
Norm = [];
FolderI = [];

for i1 = 1:size(folders,1)
    SFolder = folders{i1,1};
    
    d = dlmread([SFolder,'/',OutName]);
    
    D = [D;d];
    FolderI = [FolderI;i1*ones(size(d,1),1)];
    
end



d = D;

FInd = d(:,1);

nParam = size(d,2)-4;
al = d(:,end);

Param = d(:,2:end-3);
ParamM = Param;

% PStr = {'K_1', 'J_1', '\Gamma_1', '\Gamma_1''', 'J_2', 'J_3', 'J_z', 'D'};

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
constraint = ~isnan(Param(:,1));

% Jz being Anti-Ferro
constraint = constraint & (Param(:,7)>=0);


% constraint = constraint & (Param(:,8)>0); % 'k_1'
% cutoff = 0.05;

% % XXZ model: K = 0 and \Gamma = \Gamma'.....................
% dp = 0.05;
% constraint = constraint & (-dp<=Param(:,1)&Param(:,1)<=dp); % 'k_1'
% tt = Param(:,3) - Param(:,4); constraint = constraint & (-dp<=tt&tt<=dp); % '\Gamma_1'
% cutoff = 0.01; % 0.04

% % % Isotropic Model .........................
% dp = 0.05;
% constraint = constraint & (-dp<=Param(:,1)&Param(:,1)<=dp); % 'k_1'
% constraint = constraint & (-dp<=Param(:,3)&Param(:,3)<=dp); % '\Gamma_1'
% constraint = constraint & (-dp<=Param(:,4)&Param(:,4)<=dp); % '\Gamma'_1'
% cutoff = 0.03; % 0.04
% % %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% % %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

selCon = constraint;

chi = d(:,end-2); % Chi with Z.Z. hard constraint


chi(~selCon) = 10;

chi(chi==10) = max(chi(chi~=10));
[mchi,I0] = min(chi(1:end));
[Mchi,~] = max(chi);
chi2 = (chi-mchi)/(Mchi-mchi);

%%

% ParamM
[~,II0] = sort(chi2);
Ns = 16;
I0 = II0(1:Ns);

mParaM = min(ParamM);
MParaM = max(ParamM);
dParaM = MParaM - mParaM;
fP = find(dParaM~=0);
lp = length(fP);
% lp = size(ParamM,2);

if mod(lp,2)==0
    V = reshape(1:lp,[2,lp/2]);
else
    V = reshape([1:lp,1],[2,(lp+1)/2]);
end

IN = fP(V');
% IN = [1,2
%        2,3
%        2,4
%        2,5
%        2,6
%        8,7
% ];


le = size(IN,1);
px = ceil(sqrt(le));
py = ceil(le/px);

chit = chi2;
nchit = ((chit-min(chit))/(max(chit)-min(chit)));
selu = (nchit<cutoff);

N = 100;
[~,~,bin] = histcounts(log(chit),N);
col = jet(N+1);
ColB = col(bin+1,:);

mParam = repmat(min(Param), size(Param,1),1);
MParam = repmat(max(Param), size(Param,1),1);
NParam = (Param - mParam)./(MParam - mParam);
NParam(isnan(NParam)) = 0;

selL = find(chi2<cutoff);


FIndL = FInd(selL);
alL = al(selL);
chiL = chi2(selL);
ParamL = NParam(selL,:);

inmea = ParamL;
[coeff, score, latent, tsquared, explained] = pca(inmea);

Score = inmea*coeff;
ScoreS = NParam(I0,:)*coeff;


Ns = length(I0);

figure
plot(explained)

figure
plot3(Score(:,1),Score(:,2),Score(:,3),'bo','markerfacecolor','b')
hold on
plot3(ScoreS(:,1),ScoreS(:,2),ScoreS(:,3),'gx','linewidth',2)
text(ScoreS(:,1),ScoreS(:,2),ScoreS(:,3),num2str([1:Ns]'),'color','m','fontsize', 16,'fontweight','bold')
view(2)
rotate3d on
xlabel('P1')
ylabel('P2')
zlabel('P3')



disp('Optimial Parameter set (mean)')
mean(ParamM(selu,:))

disp('Optimial Parameter set (1:16')
mean(ParamM(II0(1:16),:))

figure

for j1 = 1:le
    subplot(px,py,j1)
    plot(ParamM(:,IN(j1,1)),ParamM(:,IN(j1,2)), 'kx'); hold on
    scatter(ParamM(selu,IN(j1,1)),ParamM(selu,IN(j1,2)),[],ColB(selu,:),'filled')
    hold on
    hold on
    text(ParamM(I0,IN(j1,1)),ParamM(I0,IN(j1,2)), num2str([1:Ns]'),'color','m','fontsize', 16,'fontweight','bold')
    xlabel(PStr{IN(j1,1)})
    ylabel(PStr{IN(j1,2)})
end

set(gcf,'position',[1          41        1600         783])


%%
% Plotting Optimization Results on Color Map.......................

LL = {[-2.,4], [-3.,2.],[-2.0,2], [-2,2],  [-2,2], [-0.5,3.5], [-1,1.5], [-1,1]};

cutoff = [inf, 0.5, 0.2,0.1, 0.07, 0.05, 0.03, 0.02, 0.01]; % 0.04


cselu = cell(length(cutoff), 1);

% selL = chi2<cutoff;
chit = chi2;
nchit = ((chit-min(chit))/(max(chit)-min(chit)));


ct = zeros(length(cutoff), size(ParamM,2));
for i1 = 1: length(cutoff)
    cselu{i1} = (nchit<cutoff(i1));
    Pt = ParamM(cselu{i1}, :);
    ct(i1, : ) = 0.1*(max(Pt) - min(Pt));
    
end

N = 100;
[~,~,bin] = histcounts(log(chit),N);
col = jet(N+1);
ColB = col(bin+1,:);

% ct = [
%     0.5, 0.5, 0.5, 0.5;
%     [0.5,0.5,0.5,0.5];
%     [0.5,0.5,0.5,0.5]/1.5;
%     [0.5,0.5,0.5,0.5]/2;
%     0.07,0.07,0.05,0.07;
%     0.07,0.07,0.05,0.07;
%     0.07,0.07,0.05,0.07];

ct(ct <=1e-2) = 1e-2;

col = jet(length(cutoff));

cRmap = cell(le, 3);

figure

for j1 = 1:le
    subplot(px,py,j1)
    
    
    I = ParamM(:,IN(j1,:));
    
    mI = min(I);
    MI = max(I);
    dI = (MI - mI)/2;
    cI = (MI + mI)/2;
    
    
    
    fac = 1.5;
    
    N = 100;
    X = linspace(cI(1)-fac*dI(1),cI(1)+fac*dI(1),N);
    Y = linspace(cI(2)-fac*dI(2),cI(2)+fac*dI(2),N);
    %     Z = linspace(cI(3)-fac*dI(3),cI(3)+fac*dI(3),N);
    
    [XX,YY] = meshgrid(X,Y);
    
    
    VV = zeros(size(XX));
    
    
    for k1 = 1:length(cutoff)
        I = ParamM(cselu{k1},IN(j1,:));
        
        [~,D] = knnsearch(I,[XX(:),YY(:)]);
        
        V = reshape(double(D<ct(k1,j1)),size(XX));
        V = imgaussfilt(V,1.8);
        
        VV = VV + V;
        
        
        %         M = contour(XX,YY,V>0.5,1,'linewidth',2);
        %         sell = M(2,:) <5;
        %         patch(M(1,sell),M(2,sell),col(k1,:))
        %         dlmwrite(sprintf('Contour_SymCo_%d_%d_%f.dat',IN(j1,:),cutoff(k1)),M)
        
        hold on
        
    end
    
    VV(VV==0) = nan;
    cRmap(j1,1:3) = {XX,YY,VV};
    
    surf(XX,YY,VV,'edgecolor','none')

    xlim([min(X), max(X)])
    ylim([min(Y), max(Y)])
    xlabel(PStr{IN(j1,1)})
    ylabel(PStr{IN(j1,2)})
    setFigStandards
    set(gca,'linewidth',2,'fontsize',16)
end



set(gcf,'position',[1          41        1920         963])


%%
% Plot Minimum
[~,II0] = sort(chi2);
I0 = II0(1);
ParamM(II0(1:16),:)
[XX,YY,Z, InsI] = readSimDataSym(folders{FolderI(I0)},FInd(I0), datafile, Range);

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









