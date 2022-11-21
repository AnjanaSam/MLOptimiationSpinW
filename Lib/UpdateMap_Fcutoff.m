% Author : Anjana M Samarakoon....................
% Date: 11/18/2022..........................
function MappingObj = UpdateMap_Fcutoff(MappingObj,nVal)

Iteration = MappingObj.Iteration;

MappingObj.AllSample = [MappingObj.AllSample;MappingObj.Positions];
MappingObj.AllFun = [MappingObj.AllFun;nVal];

disp(sprintf('Iteration %d',Iteration))
Iteration = Iteration + 1;


MappingObj.Iteration = Iteration;

