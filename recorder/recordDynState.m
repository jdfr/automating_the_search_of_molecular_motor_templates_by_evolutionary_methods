function rec = recordDynState(rec, ss, T, state)
%make sure that room is available
if rec.dynIndex>numel(rec.dyn)
  rec.dyn{end+rec.incdyn,:} = [];
  rec.dynTs(end+rec.incdyn,:) = 0;
end
%extract states and record them
index = rec.dynIndex;
rec.dynTs(index) = T;
rec.dyn{index} = state;
rec.dynIndex = index + 1;
%recordingCallback(rec, ss, T, state); %do not forget to call recordingCallback!!!

%%%%%%%%%%%%%%%
%OLD
%%%%%%%%%%%%%%%
% %make sure that room is available
% while (rec.dynIndex+numel(Ts))>numel(rec.dyn)
%   rec.dyn{end+rec.incdyn,:} = [];
%   rec.dynTs(end+rec.incdyn,:) = 0;
% end
% %extract states and record them
% indexes = rec.dynIndex+(0:(numel(Ts)-1));
% rec.dynTs(indexes) = Ts;
% rec.dyn(indexes) = mat2cell(states, size(states,1), ones(1, size(states,2)));
% rec.dynIndex = indexes(end) + 1;
% recordingCallback(rec, ss, Ts, states); %do not forget to call recordingCallback!!!
