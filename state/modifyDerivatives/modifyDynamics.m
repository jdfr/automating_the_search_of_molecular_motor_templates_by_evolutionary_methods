%adds a modification (either perturbation or enforcement of Dynamics)
%to the system.
% updateFlagVars: flag to say that affected flag vars should be set to true
function ss = modifyDynamics(ss, updateFlagVars, typemod, dynvartype, varargin)  %varargin = {indexes, startTs, endTs, modification} = modificationFields 
  num_mods = max(cellfun(@(x)size(x,1), varargin));
  for k=1:numel(ss.modificationFields)
    if (num_mods>1) && (size(varargin{k},1)==1)
      ss.(typemod).(dynvartype).(ss.modificationFields{k}) = [ss.(typemod).(dynvartype).(ss.modificationFields{k}); repmat(varargin{k}, num_mods, 1)];
    else
      if size(varargin{k},1)~=num_mods
        error('modifying dynamics for (%s).(%s): calculated number of modifications %g do not agree with number of (%s) %g', typemod, dynvartype, num_mods, ss.modificationFields{k}, size(varargin{k},1));
      end
      ss.(typemod).(dynvartype).(ss.modificationFields{k}) = [ss.(typemod).(dynvartype).(ss.modificationFields{k}); varargin{k}];
    end
  end
  if updateFlagVars
    ind   = strcmp(dynvartype, ss.dynVars);
    %varargin{1}==indexes
    ss.(ss.doubledFlagVars{ind})(varargin{1}) = true;
  end
end
