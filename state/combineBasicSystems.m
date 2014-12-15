%combine two structures complying with BasicSystems into one huge
%structure. The huge structure will inherit all specific fields from ss1,
%and will be added points, springs and modifications from ss2. New indexes
%for points and springs from ss2 are provided
function [ss1, idxPoints2, idxSprings2] = combineBasicSystems(ss1, ss2, toBeCopied)

%idxPoints1  = 1:size(ss1.pos,1);
%idxSprings1 = 1:size(ss1.springEnds,1);

idxPoints2    = size(ss1.pos,1)       +(1:size(ss2.pos,1));
idxSprings2   = size(ss1.springEnds,1)+(1:size(ss2.springEnds,1));

if isfield(ss2, 'pointIndexVars')
  reindexer = myauxReindex(idxPoints2);
  for k=1:numel(ss2.pointIndexVars)
    if ~isempty(ss2.(ss2.pointIndexVars{k}))
      ss2   = reassignField(ss2, ss2.pointIndexVars{k}, 0, reindexer);
    end
  end
end
if isfield(ss2, 'springIndexVars')
  reindexer = myauxReindex(idxSprings2);
  for k=1:numel(ss2.springIndexVars)
    if ~isempty(ss2.(ss2.springIndexVars{k}))
      ss2   = reassignField(ss2, ss2.springIndexVars{k}, 0, reindexer);
    end
  end
end

if nargin<3
  toBeCopied = [ss1.pointVars, ss1.springVars, ss1.atomVars, ss1.pShiftingVars];
end

for k=1:numel(toBeCopied)
  name = toBeCopied{k};
  ss1 = reassignField(ss1, name, 1, getFieldContents(ss2, name));
end



function fun = myauxReindex(indexes)
fun =  @(x)auxReindex(indexes, x);

% %cell matrix of indexes according to BasicSystem.dynVars = {pos,vel,m,k,r,c}; 
% indxs = {idxPoints2, idxPoints2, idxPoints2, idxSprings2, idxSprings2, idxSprings2};
% %add modifications
% for k=1:numel(ss1.modificationVars)
%   %for each modif type: perturbations & enforcements
%   modif = ss1.modificationVars{k};
%   for m=1:numel(ss1.dynVars)
%     %for each dynvar: pos vel m k r c
%     dynvar = ss1.dynVars{m};
%     %get the modification substructure
%     substruct = ss2.(modif).(dynvar);
%     %if it is not empty
%     if ~isempty(substruct.indexes)
%       %recalculate indexes
%       substruct.indexes = indxs{m}(substruct.indexes);
%       %add modification substructure
%       for n=1:numel(ss1.modificationFields)
%         fld = ss1.modificationFields{n};
%         ss1.(modif).(dynvar).(fld) = [ss1.(modif).(dynvar).(fld); substruct.(fld)];
%       end
%     end
%   end
% end