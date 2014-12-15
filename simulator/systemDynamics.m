%this function implements the dynamics of the system: given the
%system's state, it calculates its derivative
%
%ss.t MUST *NOT* BE USED!!!!! USE t INSTEAD!!!!!
function state = systemDynamics(t, state, ss)
  %let's prepare the state vars: the state is a column vector, but it
  %is convenient to rearrange it. Also, the state vector only features
  %dynamical vars, but it is convenient to rearrange them into the set
  %of all vars.
  [ss dpos] = prepareWorkingState(ss, state);
  ss.t = t;

  %point forces have two components: viscosity drag and spring forces
  dvel = calculateNaturalAcels(ss, dpos);

  %derivatives of dynamic variables. Obviously, dpos=vel and dvel=acel
  dm   = zeros(ss.numdyn_m,1);
  dk   = zeros(ss.numdyn_k,1);
  dr   = zeros(ss.numdyn_r,1);
  dc   = zeros(ss.numdyn_c,1);
  du   = zeros(ss.numdyn_u,1);
  
% %   sps = [2]; %[1 6 2];
% %   dr(sps) = 2*cos(1*ss.t)*[-1]';
%   basic = [1, 6, 61, 66, 121, 126];
%   offsets = reshape(repmat(0:9, 6, 1), 1 , [])'*pi/2;
%   n = reshape(bsxfun(@plus, basic, [0:9]'*6)', 1, []); %[2, 2+3, 8, 8+3, 14, 14+3, 20, 20+3, 26, 26+3, 32, 32+3, 38, 38+3, 44, 44+3, 50, 50+3, 56, 56+3];%[7, 7+5, 67, 67+5, 127, 127+5];%[55, 55+5, 115, 115+5, 175, 175+5]; %[67; 127];
%   f = 1;
%   f1 = -0.5;
%   f2 = f1;
%   dr(n,:) = dr(n,:)+f1*cos(f*ss.t-offsets).*(ss.t>=offsets); %bsxfun(@plus, dr(n,:), f1*cos(f*ss.t));%[f1*sin(f*ss.t) f2*cos(f*ss.t)];
  
%   dposold = dpos; dvelold = dvel; %any(any(dposold([5,6],:)-dpos([5,6],:))) || any(any(dvelold([5,6],:)-dvel([5,6],:)))
  
  %now, modify basic behaviour.
  [dpos dvel dm dk dr dc du] = nonNaturalDynamics(ss, t, dpos, dvel, dm, dk, dr, dc, du);

  %dump results on state var
  state = [dpos(:); dvel(:); dm; dk; dr; dc; du];
end




