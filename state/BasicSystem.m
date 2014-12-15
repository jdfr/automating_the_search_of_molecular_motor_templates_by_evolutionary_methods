%constructs a BasicSystem structure, which stores the simulated system's information
%
%PLEASE BE AWARE: IF THE SYSTEM NEEDS MORE STATE FIELDS (FOR EXAMPLE TO BE
%                 PASSED TO modifyXXX CALLBACKS), THE META-FIELD
%                 allStateVars ***MUST*** BE UPDATED TO ALSO INCLUDE THESE
%                 FIELDS!!!!!!!!!!!
%PLEASE BE AWARE: IS ***HIGHLY UNRECOMMENDED*** TO USE CLOSURED FUNCTIONS
%                 AS modifyXXX CALLBACKS (THAT IS TO SAY, TO USE NESTED
%                 FUNCTIONS OR ANONYMOUS FUNCTIONS REFERENCING WORKSPACE
%                 VARS): THE CLOSURED VARS WILL BE STORED IN THE CALLBACKS,
%                 WHICH ARE STORED AS STATE INFORMATION DURING RECORDING,
%                 POSSIBLY MANY TIMES
%PLEASE BE AWARE: MANY FUNCTIONS USE GENERIC CODE TO MANIPULATE DYNAMIC
%                 VARS. HOWEVER, SOME OF THEM DO NOT. SO, IF YOU WANT TO 
%                 ADD MORE DYNAMIC VARS, THE FOLLOWING FUNCTIONS MUST BE
%                 MODIFIED:
%                     -BasicSystem (add dyn var type, add fields to strcuture initialization)
%                     -addPoints (instance default params)
%                     -addSprings (instance default params)
%                     -nonNaturalDynamics (add code to handle the new var type)
%                     -prepareWorkingState (add code to handle the new var type)
%                     -makeEvolvedState (add code to handle the new var type)
%                     -systemDynamics (add code to handle the new var type)
function ss = BasicSystem(ndims, t, u, tickWidth, stepsByTick)
  if ~exist('ndims', 'var'); ndims=2; end
  if ~exist('u', 'var'); u=0; end
  if ~exist('t', 'var'); t=0; end
  if ~exist('tickWidth', 'var'); tickWidth=0.5; end
  if ~exist('stepsByTick', 'var'); stepsByTick=50; end
  %these constants are used to fold lengthy lists of assignments to struct
  %members into FOR loops
  %THE ORDERING MUST NOT BE ALTERED, THE CODE DEPENDS ON IT IN NUMEROUS WAYS 
  systemDynVars      = {'u', 't'}; %medium viscosity (1x1) AND time
  pointDynVars       = {'pos',  'vel',  'm'}; % positions (npointsXndims), velocities (npointsXndims), masses (npointsX1)
  springDynVars      = {'k',  'r',  'c'}; %K constants (nspringsX1), rest lengths (nspringsX1), friction constants (nspringsX1)
  dynVars            = {pointDynVars{:}, springDynVars{:}};%, systemDynVars{1}};
  springVars         = {'springEnds', springDynVars{:}}; 
  springDefaultVals  = {zeros(1,2), 1, 1, 0}; %default values for members of springVars
  springIndexVars    = {}; %vars whose content are spring indexes
  pointVars          = {pointDynVars{:}};
  pointDefaultVals   = {zeros(1,ndims), zeros(1,ndims), 1}; %default values for members of pointVars
  pointIndexVars     = {springVars{1}}; %vars whose content are point indexes
  allStateVars       = {systemDynVars{:},    ...
                        pointVars{:},        springVars{:}, ...
                        'u', 'pointDefaultVals',  'springDefaultVals', 'dynamic', ...
                            'brownianFractor', 'useAnisodamp', 'useSAP' };
  fieldSets             = {'pointDynVars', 'springDynVars', ...
                           'dynVars', ...
                           'pointVars', 'springVars', 'pointIndexVars', 'springIndexVars', ...
                           'allStateVars', 'fieldSets'};
  %helper initial values
  complexReshaper = @(strcells,args) reshape([strcells;args],[],1);
  %simpleReshaper  = @(strcells) complexReshaper(strcells, cell(1,numel(strcells)));
  %properly assign initial dimensions to vars
  dynvars       = complexReshaper(dynVars, {zeros(0,ndims), zeros(0,ndims), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1)});%, u});
  dynFlags      = complexReshaper(dynVars, {true, true, false, false, true, false});
  
  %create the STRUCTURE
  ss = struct(...
    ... %STATE VARS
    'tick', tickWidth, 'stepsByTick', stepsByTick, ...
    't', t, 'u', u, ...                   %time, u, 
    'springEnds', zeros(0,2), ... %array of point indexes for springs' ends: nspringsX2
    dynvars{:}, ...
    'dynamic', struct(dynFlags{:}), ...
    ... %OPTIONS VARS
    'chgStep', [], ...
    'postProcess', [], ...
    'brownianFractor', 0, ...
    'useAnisodamp', false, ...
    'useSAP', false, ...
    ... %The two following variables are structs of structs. Each is a struct
    ... %containing a struct for each dynamical variable. Each modification is
    ... %related not to the variable, but its derivative
    ... %these constants are used to fold lengthy lists of assignments to struct
    ... %members into FOR loops
    ... %THE ORDERING MUST NOT BE ALTERED, THE CODE DEPENDS ON IT IN NMEROUS WAYS 
    'dynVars',                {dynVars}, ...
    'pointVars',              {pointVars}, ...
    'pointDefaultVals',       {pointDefaultVals}, ...
    'springVars',             {springVars}, ...
    'springDefaultVals',      {springDefaultVals}, ...
    'pointDynVars',           {pointDynVars}, ...
    'springDynVars',          {springDynVars}, ...
    'pointIndexVars',         {pointIndexVars}, ...
    'springIndexVars',        {springIndexVars}, ...
    'allStateVars',           {allStateVars}, ...
    'fieldSets',              {fieldSets}, ...
    ... THESE FIELDS ARE FOR INTERNAL USE; THEY ARE UPDATED BY prepareSimulationStructure
    'npoints',  [], ... %number of points
    'nsprings', [], ... %number of springs
    'ndims',    ndims, ... %number of dimensions
    ... %working var:
    ... %     (i,j) pos == 1 => i-th node is first  end of j-th spring
    ... %     (i,j) pos ==-1 => i-th node is second end of j-th spring
    'springsMatrix', [] ...
  );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rdynfields      = cellfun(@(x){'rdyn', x}, {'radapt', 'r0',         'rtb',  'rt', 'rspan'}, 'uniformoutput', false);
rdynfieldsvals  =                          {false,    'PREVIOUS_r', 0,       1,   0};
newspringfields        = rdynfields;%     rgenfields];% newspringfieldsk];
newspringdefvals       = rdynfieldsvals;% rgenfieldsvals];% newspringdefvalsk];
for k=1:numel(newspringfields)
  if (~ischar(newspringdefvals{k})) || (~strmatch('PREVIOUS_', newspringdefvals{k}))
    res = repmat(newspringdefvals{k}, size(ss.k), 1);
  else
    res = ss.(newspringdefvals{k}(10:end));
  end
  ss   = reassignField(ss, newspringfields{k}, 0, res);
end

ss.rdyn.rchg = false;

ss.springVars         = [ss.springVars,         newspringfields];
ss.springDefaultVals  = [ss.springDefaultVals,  newspringdefvals];

newpointvars    = cellfun(@(x){'stick', x}, {'rad', 'allr', 'penk'}, 'uniformoutput', false);
newpointdefvals =                           {0,     0,      1};

for k=1:numel(newpointvars)
  ss   = reassignField(ss, newpointvars{k}, 0, repmat(newpointdefvals{k}, size(ss.m)));
end

ss.pointVars          = [ss.pointVars,        newpointvars];
ss.pointDefaultVals   = [ss.pointDefaultVals, newpointdefvals];

ss.allStateVars       = [ss.allStateVars,     {'rdyn', 'stick'}];


end
