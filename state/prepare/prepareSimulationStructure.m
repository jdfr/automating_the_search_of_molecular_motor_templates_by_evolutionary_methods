%makes a helper structure for the simulation. It precalculates all derived
%fields. It simply calls other more specialized functions
function ss = prepareSimulationStructure(ss)
%to be used before simulating if the number of dimensions, points or
%springs changes 
ss= makeSimulationAmounts(ss); 
%to be used before simulating if the the number of points or spring
%changes, or the connectivity of any spring changes
ss =  makeSpringEndsMatrix(ss);

if ss.useSAP && ( (~isfield(ss, 'sap')) || ((numel(ss.sap.idxsx)/2)~=size(ss.pos,1)) )
  if isfield(ss, 'updateSAP')
    ss.sap = ss.updateSAP(ss);
  else
    ss.sap = createSweepAndPrune(ss.pos, ss.stick.allr);
  end
end