function executePicassoTask(nomdir, inputFileNameTemplate, outputFileNameTemplate, numTask, mypath)

%#function molMot3DConf
%#function pivotFunMolMot3DConf
%#function orientFunMolMot3DConf
%#function evaluarMolMots3D
%#function createRandomFolds
%#function mutateFold
%#function ANMdecomposition
%#function collisionsWithRow
%#function convertCompressedPos
%#function evaluarMolMot
%#function doEvaluationGaGeneric
%#function manipulateGenome2D
%#function manipulateGenome3D
%#function prepareLegsSimulation
%#function rotateLegs
%#function snapToGrid
%#function structure2walker
%#function transformCorrectParams
%#function validChain
%#function walkerMachine
%#function DAConvertGenome
%#function DADevelopGenome
%#function DAProcessGenes
%#function DAcreateRandomGenome
%#function DAgenomesAreEqual
%#function DAmakeMutationConfiguration
%#function DAnumCellsInGenome
%#function DAsplitGenome
%#function DynamicAtomSystem
%#function addAtomComplex
%#function calculateOverload
%#function clearOldPerturbations
%#function createAtomSprings4P1S
%#function decodeBest
%#function magicallyMovedPoints
%#function perturbateR
%#function postProcessPoph
%#function removeAtoms
%#function sowSeed
%#function MemRecorder
%#function makeNewRecorder
%#function recordAllState
%#function recordDynState
%#function recordingCallback
%#function calculateAcels
%#function calculateEnergy
%#function calculatePotentialEnergy
%#function prepareWorkingState
%#function simulateSystem
%#function simulateSystemIntercalated
%#function simulateSystemIntercalatedCoarse
%#function simulateSystemUntilStability
%#function createSweepAndPrune
%#function doSweepAndPrune
%#function BasicSystem
%#function combineBasicSystems
%#function copyFields
%#function findSpringIndexes
%#function findSpringsForPoints
%#function getBasicSystemFieldSets
%#function getFieldContents
%#function reassignField
%#function addPoints
%#function addSprings
%#function auxAssignment
%#function auxReindex
%#function auxRemoveModifications
%#function copyPoints
%#function copySprings
%#function findIndexesAfterRemove
%#function removePoints
%#function removeSprings
%#function removeUnconnectedPoints
%#function dumpEvolvedState
%#function makeDynamicVars
%#function makeEvolvedState
%#function makeSimulationAmounts
%#function makeSpringEndsMatrix
%#function prepareSimulationStructure
%#function alwaysTrue
%#function capPrecision
%#function pointDistancesAtTime
%#function calculateRanges
%#function normalDist
%#function pickOption
%#function affinTransform
%#function closestCircleLineIntersection
%#function curvature
%#function distanceMatrix
%#function distanceMatrixSquared
%#function drawAsImage
%#function intersectionOfLines
%#function longSegments
%#function pdist1D
%#function pointPolyDist
%#function imcapture
%#function mprintf
%#function stoploop
%#function vprintf
%#function mergesorted
%#function filenamedirrec
%#function my_filenamedirrec
%#function any2str
%#function array2cell
%#function assignSubStruct
%#function assignSubStructOld
%#function getSubStruct
%#function getSubStructOld
%#function makeEmptyStructure
%#function makeListStructure
%#function repmatStruct
%#function showError


numTask = str2double(numTask);

% % fprintf('Mira0: %s\n', pwd);
% 
% cd ../..
% 
% % fprintf('Mira1: tralari tralara\n');
% % fprintf('Mira2: %s\n', pwd);
% % fprintf('Mira3: %s\n', nomdir);
% % fprintf('Mira4: %s\n', inputFileNameTemplate);
% % fprintf('Mira5: %s\n', outputFileNameTemplate);
% % fprintf('Mira6: %d\n', numTask);
% % exit

inputFile  = [nomdir sprintf(inputFileNameTemplate, numTask)];
outputFile = [nomdir sprintf(outputFileNameTemplate, numTask)];

if exist('mypath', 'var') && (~isempty(mypath))
  addpath(mypath);
end

load(inputFile);

%this is in order to not copy evparams in multiple separate files
if ischar(task.InputArguments{1}) %#ok<NODEF>
  load([nomdir task.InputArguments{1}])
  task.InputArguments{1} = evparams;
end


output = cell(task.NumberOfOutputArguments, 1);

ok = true;

try
  [output{:}] = task.Function(task.InputArguments{:}); %#ok<NASGU>
catch ME
  err = showError(ME); %#ok<NASGU>
  ok = false;
end

if ok
  save(outputFile, 'output', '-v7');
else
  save(outputFile, 'err', 'task', '-v7');
end

exit


% function doAddPath(nomdir)
% srcdir = nomdir;
% seps   = find(srcdir==filesep);
% srcdir = srcdir(1:seps(end-1));
% load([srcdir 'estado.mat']);
% 
% jobArgs = params.terclus.jobArgs;
% 
% for k=1:numel(jobArgs)
%   if strcmpi(jobArgs{k}, 'PathDependencies')
%     pth  = jobArgs{k+1};
%     pth  = [pth, repmat({pathsep}, size(pth))]; %#ok<AGROW>
%     pth  = reshape(pth', 1, []);
%     pth  = horzcat(pth{:});
%     addpath(pth);
%     return
%   end
% end
% 
% error('This should not happen!!!!');
% 
% 
