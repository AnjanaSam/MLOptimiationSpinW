% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................

% Note: Before running this script, make sure that you have correctly setup
% the "calculate_sw2_RealS02_KNiAsO4_Opt04.m"

% Experimental Data............
datafile = './kniaso4_13mev_300hz_5K_Bkg.bin';
% Project Folder. Define a unique name for each time you run this script.
SFolder = './DirectRunSym_kniaso4_1bK';
% Project File. This contains all the information about the optimiation
% process. Define a unique name for this as well.
outSwam = 'DirectRun_RealS_Sym_kniaso4_1bK.mat';
% cutoff sweep parameters
PcutoffMax = 1; % Max
PcutoffMin = 0.001; % Min
nL = [30,30]; % # Steps [a,b] a - # steps while lowering cut-off from Max to Min, b - # steps at Min
% Region of interest in Energy Transfer (meV)
Range = [0.5,11]; 
% # of spherical samples for powder average
NQ = 300; 
% Smoothing co-efficient. (Set it to 0)
smf = 0.0; 


% Number of Particles - Number of simulations per iteration
% This should be a multiple of number of logical cores of the paralell
% session.Ususally it is 12 in analysis.sns.gov. Might be different for
% some other machine.
numParticles = 24; 

% Boundary of Parameters: lb - lower bound, ub - upper bound 
% Order: K, J, G, Gp, J2, J3, Jz, A0
lb = [0,  -2.5,    0,     0,  -1.5, -1.5, 0,  -1.5];
ub = [0,   1.5,    0,     0,   1.5,  3.5,  1,  1.5];

% Weights for cost-function components
% e.g: chi^2 = n1*chi^2_{Q1} + n2*chi^2_{Q1}
% ID = [n1, n2, ...]
ID = [1];




%........Script Starts here.........................
%........Do Not Chnage Below this line..............
%...................................................

OutName= 'log_RealSMag.dat'; % cOut{i1} = [chi, chiOri, AAlpha];

% Setting Radom Seed................
c = clock;
seed = c(6);
rng(seed);


% Set Path to library.....................
% folder = fileparts(which('../spinw3_R1343/spinw3_R1343/'));
% Add that folder plus all subfolders to the path.
addpath('../spinw3_R1343/spinw3_R1343',...
    '../spinw3_R1343/spinw3_R1343/dat_files',...
    '../spinw3_R1343/spinw3_R1343/external',...
    '../spinw3_R1343/spinw3_R1343/external/chol_omp',...
    '../spinw3_R1343/spinw3_R1343/external/colormap',...
    '../spinw3_R1343/spinw3_R1343/external/eig_omp',...
    '../spinw3_R1343/spinw3_R1343/external/fminsearchbnd',...
    '../spinw3_R1343/spinw3_R1343/external/mtimesx',...
    '../spinw3_R1343/spinw3_R1343/external/x3d_version3g',...
    '../spinw3_R1343/spinw3_R1343/external/x3d_version3g/functions',...
    '../spinw3_R1343/spinw3_R1343/external/x3d_version3g/functions_xml',...
    '../spinw3_R1343/spinw3_R1343/swfiles',...
    '../spinw3_R1343/spinw3_R1343/swfiles/usefull',...
    './Lib/');

nPara = length(ub);




% Region_Name = 'SelectedRegion3.dat';

ND = length(ub);

% Parical Mapping Parameters......................
% Factor to adjust number of function evaluations
fac = 2^ND; % 10

% Anneal Parameters............................
nItr = 10; % 10  iteraction per level
% nL = [20,20]; % [40,20]; % 10 linear levels 5 constant levels

% Selected Region
% 0: Selected   1: Avoided
SReg = [(lb+ub)/2,0];
% SReg = dlmread(Region_Name);


% Restate Parameter...............................
restart = 1;

if ~restart
    
    Plevel = [linspace(PcutoffMax,PcutoffMin,nL(1)),PcutoffMin*ones(1,nL(2))];
    
    
    
    if ~exist(SFolder)
        mkdir(SFolder);
    end
    if ~exist([SFolder,'/',OutName])
        fid = fopen([SFolder,'/',OutName], 'wt' );
        fclose(fid)
    end
    
    
    if exist(outSwam)
        mf  = matfile(outSwam);
        MappingObj = mf.MappingObj;
        MappingObj.Plevel = [MappingObj.Plevel,Plevel];
    else
        MappingObj = InitializeParticleSwarm_FcutoffSR(ub, lb, nPara, numParticles,PcutoffMax,Plevel,fac, SReg);%*******************<<<<<<<<<<<<<<
        % Load Data from Previous Files.......
        folders = {
%             './DirectRunSym_kniaso4_0bS';
%                 './DirectRunSym_kniaso4_0bK';
            };
        
        AP = [];
        AF = [];
        
        for z1 = 1:length(folders)
            SFoldert = folders{z1};
            d0 = dlmread([SFoldert,'/', OutName]);
            
            AP = [AP; d0(:,2:end-3)];
            AF = [AF; d0(:,end-2)];
        end
        
        MappingObj.AllSample = AP;
        MappingObj.AllFun = AF;
        
        
        %..........................................
    end
    save(outSwam,'MappingObj')
    stI1 = 0;
    stI2 = 0;
else
    mf  = matfile(outSwam);
    MappingObj = mf.MappingObj;
    Plevel = MappingObj.Plevel;
    stI1 = floor(MappingObj.Iteration/nItr);
    stI2 = mod(MappingObj.Iteration,nItr);
end

save(outSwam,'MappingObj')

try
    for i1 = stI1+1:length(Plevel)
        MappingObj.Pcutoff = Plevel(i1);
        for j1 = stI2+1:nItr
            
            MappingObj = GetNewPosition_FcutoffKNNS_Constrain04(MappingObj, [], [], ID);%*******************<<<<<<<<<<<<<<
            pos = MappingObj.Positions;
            %             nVal = getminChi2M(pos,SFolder);
            nVal = getChi2M_RealS_Sym_KNiAsO4R04(pos,SFolder,1,Range, NQ, datafile, smf, OutName);
            
            MappingObj = UpdateMap_Fcutoff(MappingObj,nVal);%*******************<<<<<<<<<<<<<<
        end
        MappingObj = Update_ParamLimits(MappingObj,ID,5);
        save(outSwam,'MappingObj')
    end
    
catch ME
    save(outSwam,'MappingObj')
    
    rethrow(ME)
end

save(outSwam,'MappingObj')
















