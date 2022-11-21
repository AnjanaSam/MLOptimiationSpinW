
function [out,magStrF] = calculate_sw2_RealS02_KNiAsO4_Opt04(Param,Ind,SFolder,q,de, NQ, EvalFlag)

if nargin < 7
    EvalFlag = 1;
end


%.................................................
% Creating Spinwave Model.........................
%..................................................
KNiAsO4 = spinw('./As1 K1 Ni1 O4_Ni_Only.cif');
KNiAsO4.unit_cell.S(1) = 1;
% KNiAsO4.unit_cell.r = [1/3, 2/3, 0]';
% plot(KNiAsO4)


KNiAsO4.gencoupling('maxDistance',12)
% KNiAsO4.table('bond',[])




% This is Correct....................
NN = [-0.4084,-0.7071,-0.5773;-0.4084,0.7071,-0.5773;0.8164,0,-0.5775];

% getTransMat(NN(:,1))*Rz(phi)*getTransMat(V2)'
% j1x = getTransMat(NN(:,1))'*diag([1,0,0]*J1)*getTransMat(NN(:,1));



% j0x = diag([1,0,0]*J1);
% j0y = diag([0,1,0]*J1);
% j0z = diag([0,0,1]*J1);



K = Param(1); % 1; 
J = Param(2); % 0.5;
G = Param(3); % 0;
Gp = Param(4); % 0.00;

J2 =Param(5); % -0.6;
J3 = Param(6); % -0.6;
Jz = Param(7); % 0.1;
A0 = Param(8);




%J1=-1.0*abs(0.54/(1+2*alpha+beta)); % BASED ON CW temperature
% J1=0.54*1.625/(1+2*alpha+beta);
% j0x = [  J, Gp,  G;
%         Gp, K+J, Gp;
%          G, Gp,   J];
j0x = [  K+J, Gp,  Gp;
        Gp, J, G;
         Gp, G,   J];
% j0x = [  J, G,  Gp;
%         G, J, Gp;
%          Gp, Gp,   K+J];     


j1x = NN'*j0x*NN;
% j1y = NN'*j0y*NN;
% j1z = NN'*j0z*NN;


%%%%%%%%%%%%%%%%
KNiAsO4.addmatrix('label','J1','value',j1x,'color','red')     % 1.2
KNiAsO4.addmatrix('label','J2','value',J2,'color','blue') % 0.1
KNiAsO4.addmatrix('label','J3','value',J3,'color','green') % -1.3
KNiAsO4.addmatrix('label','Jz','value',Jz,'color','yellow')
KNiAsO4.addmatrix('label','A','value',diag([0 0 -A0]))  % easy-axis

KNiAsO4.addcoupling('mat','J1','bond',1)
KNiAsO4.addcoupling('mat','J2','bond',2)
KNiAsO4.addcoupling('mat','J3','bond',3)
KNiAsO4.addcoupling('mat','Jz','bond',7)
KNiAsO4.addaniso('A')

%defining the magnetic structure 
Sseq = [-1,1, 1, -1, 1, -1]; % ZZ
Kin = [0.5 0.5 0];
S0 = [1,1,0]'*Sseq;
KNiAsO4.genmagstr('mode','helical','k',Kin,'n',[0 0 1], 'S',S0,'nExt',[2 2 1]);

if ~EvalFlag
    plot(KNiAsO4,'range',[2 2 1]); hold on;
end

% display('Ground state energy (meV/spin)')
E1 = KNiAsO4.energy;
% plot(Na2Ni2,'range',[2 2 1]); hold on;

% Optimize for Magntic Structure......................
KNiAsO4.optmagsteep('nRun',1e6);
    E2 = KNiAsO4.energy;
%     E2-E1
%     plot(Na2Ni2,'range',[2 2 1]); hold on;

if ~EvalFlag
    plot(KNiAsO4,'range',[2 2 1]); hold on;
end
swpref.setpref('usemex',true)

if EvalFlag
    AFsqpowspec = KNiAsO4.powspec(q,'Evect',[0,de],'nRand',NQ,'formfact',true,'hermit',false,'imagChk',false);
    % AFsqpowspec = sw_instrument(AFsqpowspec,'dE',0.5,'dQ',0.02,'Ei',15000,'thetaMin',0);
    
    out = AFsqpowspec.swConv;
    
    St = KNiAsO4.magstr.S;
    
    % Condition for Z.Z. Ordering.....
    % Might not True for another magnetic structure....
    tt = diff(unique(roundn(St,-3)','rows'),[],1);
    magStrF = (size(tt,1)==1)&&round(norm(tt))==2;
    
    writeBin([SFolder,'/Model_',num2str(Ind),'.bin'],{q,de,out, E1,E2,St})
    
end
% writeBin([SFolder,'/Model_',num2str(Ind),'.bin'],{})


% % sw_plotspec(AFsqpowspec,'mode',3,'dE',0.3); hold on; axis([0.2 2.8 0 13]); caxis([0 1]); shading interp; colormap(jet(256));
% % 
% % AFsqpowspec.swConv=AFsqpowspec.swConv.*Instr;
% % 
% % figure(101); sw_plotspec(AFsqpowspec,'mode',3,'dE',0.3); axis([0.2 2.75 1 14]); caxis([0 1]); shading interp; colormap(jet(256));