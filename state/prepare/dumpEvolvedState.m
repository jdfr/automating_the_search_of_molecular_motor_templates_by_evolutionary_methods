%dump 'evolvedState' into the structure's dynamical vars
function ss = dumpEvolvedState(ss, evolvedState)
  fn = fieldnames(evolvedState);
  for k=1:numel(fn)
    ss.(fn{k}) = evolvedState.(fn{k});
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %common components unfolded, to speed the code up
%   dims = size(ss.pos,2);
%   ss.pos(ss.dynamical_p,:) = reshape(evolvedState(ss.selector.pos(1,1):ss.selector.pos(end,2)), [], dims);
%   ss.vel(ss.dynamical_p,:) = reshape(evolvedState(ss.selector.vel(1,1):ss.selector.vel(end,2)), [], dims);
%   if any(ss.dynamical_m)
%     ss.m(ss.dynamical_m,:)   = evolvedState(ss.selector.m(1,1):ss.selector.m(end,2));
%   end
%   if any(ss.dynamical_k)
%     ss.k(ss.dynamical_k,:)   = evolvedState(ss.selector.k(1,1):ss.selector.k(end,2));
%   end
%   if any(ss.dynamical_r)
%     ss.r(ss.dynamical_r,:)   = evolvedState(ss.selector.r(1,1):ss.selector.r(end,2));
%   end
%   if any(ss.dynamical_c)
%     ss.c(ss.dynamical_c,:)   = evolvedState(ss.selector.c(1,1):ss.selector.c(end,2));
%   end
% %   if ss.dynamical_u
% %     ss.u(ss.dynamical_u) = evolvedState(ss.selector.u(1));
% %   end
%   %if the implementation uses more dynamical vars, they are processed the
%   %old way
%   if (numel(ss.dynVars)>7)
%     for k=8:numel(ss.dynVars)
%       campo      = ss.dynVars{k};
%       dynamicals = ss.doubledFlagVars{k};
%       dims       = size(ss.(campo),2);
%       ss.(campo)(ss.(dynamicals),:) = reshape(evolvedState(ss.selector.(campo)(1,1):ss.selector.(campo)(end,2)), [], dims);
%     end
%   end
% 
% %   %old, fully automatized way
% %   for k=1:numel(ss.dynVars)
% %     campo      = ss.dynVars{k};
% %     dynamicals = ss.doubledFlagVars{k};
% %     dims       = size(ss.(campo),2);
% %     ss.(campo)(ss.(dynamicals),:) = reshape(evolvedState(ss.selector.(campo)(1,1):ss.selector.(campo)(end,2)), [], dims);
% %   end
% end
