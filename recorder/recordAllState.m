%record all system state
function rec = recordAllState(rec, ss)
%if ss.t is the same as the last rec.notdynTs, overprint it
if (rec.notdynIndex>1) && (rec.notdynTs(rec.notdynIndex-1)==ss.t)
  rec.notdyn{rec.notdynIndex-1} = extractAllSystemState(rec, ss);
  return
end
%make sure that room is available
if rec.notdynIndex>numel(rec.notdyn)
  rec.notdyn{end+rec.incnotdyn,:} = [];
  rec.notdynTs(end+rec.incnotdyn,:) = 0;
end
%extract state and record it
rec.notdynTs(rec.notdynIndex) = ss.t;
rec.notdyn{rec.notdynIndex} = extractAllSystemState(rec, ss);
rec.notdynIndex = rec.notdynIndex+1;
%recordingCallback(rec, ss); %do not forget to call recordingCallback!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%extracts info from BasicSystem into a structure
function st = extractAllSystemState(rec, ss)
if isempty(rec.emptyState)
%   %THIS IS VALID ONLY IF NO SYSTEM STATE IS ALLOWED TO BE NESTED IN
%   %SUBSTRUCTURES
%   aux = reshape([ss.allStateVars;cell(1,numel(ss.allStateVars))],[],1);
%   rec.emptyState = struct(aux{:});
  rec.emptyState = struct;
  for k=1:numel(ss.allStateVars)
    rec.emptyState = reassignField(rec.emptyState, ss.allStateVars{k}, 0, []);
  end
end
% st = rec.emptyState;
% for k=1:numel(ss.allStateVars)
%   %THIS IS VALID ONLY IF NO SYSTEM STATE IS ALLOWED TO BE NESTED IN
%   %SUBSTRUCTURES
%   campo = ss.allStateVars{k};
%   st.(campo) = ss.(campo);
% end
st = copyFields(ss, rec.emptyState,  ss.allStateVars);
if isfield(ss, 'sap')
  st.sap = ss.sap;
end


