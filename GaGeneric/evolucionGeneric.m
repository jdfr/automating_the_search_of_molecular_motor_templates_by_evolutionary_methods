%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G indexNewG terclus P parentdata] = evolucionGeneric(params,P,generacion,rangeid,generationInRange,terclus,mejorEnConserva)
%
%  range: 2-tuple marking the range which will be replaced. The new
%  individuals do not immediately replace the population in that range, but
%  are stored in G.
%  The indexes of truly new individuals are in indexNewG (the others are
%  just straightforward copies of individuals in P, i.e., they haven't
%  undergone any mutation at all). IMPORTANT: indexNewG is relative to G,
%  not P!!!!
%  Their values 'b' and 'ra' for new individuals are stored in 'bG' and
%  'raG' respectively: the new individuals have initially 0 (as they must
%  be evaluated), but old individuals' values are straightforwardly copied.
%
%  All this stuff is a nuissance, but it enables the code to overlap
%  execution for several subpopulations, thus being able to use the cluster
%  for a subpopulation while applying genetic operators to another
%  subpopulation
%
%   prob : probabilidades de cada operador
%
%      1   elitismo             =  
%      2   op_alteracion        A simbolo  posicion  
%      3   op_dupaleatoria      R segmento posicion
%      4   op_dupnivel          L segmento posicion
%      5   op_dupsecuencia      C segmento posicion
%      6   op_eliminacion       D posicion
%      7   op_insercion         I simbolo  posicion
%
% el patrón utilizado en el fichero arbol.txt es:
% generacion progenitor operador_segmento_posicion
% el patrón operador_segmento_posicion podrá aparecer repetido tantas veces
% como operadores se hayan aplicado sobre el individuo.
% Ejemplo un operador: 1 186 N[FGF-GGF]3
% Ejemplo varios operadores: 1 16 C+8N[+GFF++F]3
% Cuando un operador se aplica al individuo pero no modifica la cadena que
% lo representa se representará por el símbolo '=' y el parámetro posición no se almacena. 
% Ejemplo: 1 95 S=I+8 (generacion progenitor operador_=_operador_segmento_posicion)
% Por último, el patrón para el operador elitista es:
% generacion progenitor =
% Ejemplo: 1 207 =

elitista                  = (rand<=params.elitism);       % estrategia elitista (mantener mejor)
elitista                  = elitista && (terclus.ranges(rangeid,1)==1); %only if the range starts in the 1st place
  
range                   = terclus.ranges(rangeid,:);

rangesExpanded          = zeros(size(P.genome));
for k=1:size(terclus.ranges,1)
  rangesExpanded(terclus.ranges(k,1):terclus.ranges(k,2)) = k;
end

writeRange = isfield(params, 'writeRangeID') && params.writeRangeID;
if ~isfield(params, 'jiggleSelection')
  params.jiggleSelection = 0;
end

rewindNans   = params.rewindNans;
nansReplaced = params.nansReplaced;

if rewindNans
  %if individuals with fitness==nan are to be rewound, check whether there
  %are Nan fitnesses. If so, revert the individuals to the previous version
  Preverted = params.Preverted;
  rwdL = islogical(Preverted);
  if (rwdL && any(Preverted)) || ((~rwdL) && (~isempty(Preverted)))
    P  = assignSubStruct(P, Preverted, 'extensive', params.Pold);
  end
end

if nansReplaced
  numnans = sum(isnan(P.fitness));
end

numG = range(2)-range(1)+1;
useQuota = (~strcmp(params.newsQuotaMode, 'none')) || (nansReplaced && (numnans>0));
if useQuota
  %If there is a quota of new individuals, get the amount.
  switch params.newsQuotaMode
%     case 'none'
%       useQuota = false;
    case 'plain'
      quota = min(numG, round(params.newsQuota*numG));
    case 'rand'
      minq  = max(0, min(numG, round(min(params.newsQuota)*numG)));
      maxq  = max(0, min(numG, round(max(params.newsQuota)*numG)));
      quota = rand*(maxq-minq)+minq;
    case 'custom'
      quota = params.newsQuota(generacion, rangeid, generationInRange, range, P, params);
      quota = max(0, min(numG, quota));
    case 'none'
      quota = 0;
    otherwise
        error('Parameter params.newsQuotaMode==%s not understood!!!', any2str(params.newsQuotaMode));
  end
  if nansReplaced
    %this can operate on all individuals with NaN fitness, or only 
    %in those having no sane ancestors (in the latter case, the ones 
    %being the product of a bad mutation would be filtered in by the
    %rewindNans process)
    quota = min(numG, quota+numnans); 
  end
  numG  = numG-quota;
end

if numG>0
  switch params.selection.mode
    case 'PACO'
      %classical way: stocastically select individuals according to their
      %fitnes
      rango                   = calculateRangoSeleccion(P,P.fitness,params.presion,params);
      idxInP                  = seleccion(rango,numG,params);
    case 'FUSS'
      idxInP                  = FUSSfun(P, P.fitness, params, numG, []);
    otherwise
      error('params.selection.mode==%s not understood!!!\n', any2str(params.selection.mode));
  end
else
  idxInP = [];
end

parentdata  = struct;

if useQuota
  %if there is a quota, make the necessary arrangements to deal with it
  parentdata.notNew   = [false(quota, 1); true(numG, 1)];
  G                   = repmatStruct(params.indivStruct, [quota,1]);
  G.genome(1:quota)   = params.generarGenomas(quota, params);
  if numG>0
    G                 = assignSubStruct(G, quota+[1 numG], 'intensive', P, idxInP, 'extensive');
  end
  idxInP              = [zeros(quota, 1); idxInP(:)];
else
  G           = getSubStruct(P, idxInP, 'extensive');
end

% G           = repmatStruct(params.indivStruct, size(P.genome));
% idxInP      = zeros(size(G.genome));
maskNew     = false(size(G.genome));
if (elitista)
  G.genome{1} = mejorEnConserva.data.genome{1};
  maskNew(1)  = false; %true; %the best organism must revalidate itself generation after generation, since different random seeds might give him different fitnesses
  idxInP(1)   = mejorEnConserva.num;
  if useQuota
    %make sure that elitist individuals are not marked as newly generated
    parentdata.notNew(1) = true;
  end
  if writeRange
    params.writeGenealogy(generacion,rangeid,1,'=',mejorEnConserva.gen, mejorEnConserva.rind,mejorEnConserva.num);
  else
    params.writeGenealogy(generacion,        1,'=',mejorEnConserva.gen,                      mejorEnConserva.num);
  end
end

if rewindNans
  parentdata.idxInP = idxInP;
end

lggng = params.logInMutation;%params.logInEval;
if lggng
  zz = [params.absdir 'zmutations.txt'];
  [f msg]=fopen(zz, 'w');
  if ~isempty(msg)
    error('We have not been able to open <%s>\nCause: <%s>\n', zz, msg);
  end
  fprintf(f, 'NUMBER OF GENOMES TO "MUTATE": %d\n', length(G.genome));
  fclose(f);
end

for j=elitista+1:length(G.genome)      % crear individuos a partir de la anterior
    i    = idxInP(j);

    if lggng; f=fopen(zz, 'a'); fprintf(f, 'GENOME %03d (%03d in prev. pop.), BEFORE AND AFTER:\n%s\n', j, i, params.genome2str(G.genome{j})); fclose(f); end
    
    if useQuota && (i==0)
      %if the indivudals has been generated afresh, it has no parent
      maskNew(j)  = true;
      cad         = '';
      ranOld      = 0;
      genOld      = 0;
    else
      [G.genome{j} cad maskNew(j)] = params.mutateGenome(j, G, params, P, idxInP, maskNew);
      %[G.genome{j} cad maskNew(j)] = mutateGenomeSprings(G.genome{j}, params, numMutations(j), maxLongsImpl);

      if lggng; f=fopen(zz, 'a'); fprintf(f, '%s\nCHANGE: %s\n', params.genome2str(G.genome{j}), cad); fclose(f); end
    
      ranOld = rangesExpanded(i);
      genOld = generationInRange(ranOld);
      
      if rewindNans && ( (rwdL && Preverted(i)) || ((~rwdL) && any(Preverted==i)) )
        %make sure that the condition of reverted individuals is annotated
        cad = ['[REVERTED]' cad]; %#ok<AGROW>
      end
    end
    indexNewInd = j+range(1)-1;
      
    if ~isempty(cad)
        %añadir al fichero de texto: generación, indice del nuevo
        %individuo, indice del viejo individuo, cambios
        if writeRange
          params.writeGenealogy(generacion,rangeid,indexNewInd,cad,genOld,ranOld,i);
        else
          params.writeGenealogy(generacion,        indexNewInd,cad,genOld,       i);
        end
    else
        %añadir al fichero de texto
        if writeRange
          params.writeGenealogy(generacion,rangeid,indexNewInd,'=',genOld,ranOld,i);
        else
          params.writeGenealogy(generacion,        indexNewInd,'=',genOld,       i);
        end
    end
end

if lggng; delete(zz); end; %f=fopen(zz, 'a'); fprintf(f, 'about to delete this file\n'); fclose(f); delete(zz); end

if params.evaluateAlways && (~isfield(G, 'rasterDev'))
  %all of them must be re-evaluated, since a good fitness might be a fluke
  %from a good configuration from a random perturbation
  maskNew(1:end) = true; 
  indexNewG      = 1:numel(maskNew);
else
  if isfield(G, 'raster')
    %ARREGLAR ESTO!!!!!
    %maskNew         = maskNew | (cellfun('prodofsize', G.raster)==0);
  elseif isfield(G, 'rasterDev')
    maskNew         = maskNew | (cellfun('prodofsize', G.rasterDev)==0);
  end
  %get the indexes for new individuals
  indexNewG       = find(maskNew);
end


%put default values in all new individuals
fieldsNoGenome  = fieldnames(G);
fieldsNoGenome(strcmp('genome', fieldsNoGenome)) = []; %do not assign the genome field
% G = assignSubStruct(G, ~maskNew, 'extensive', P, idxInP(~maskNew), 'extensive', fieldsNoGenome);
G               = assignSubStruct(G, maskNew, 'extensive', params.indivStruct, 1, 'extensive', fieldsNoGenome);

if elitista
  G      = assignSubStruct(G, 1, 'extensive', mejorEnConserva.data, 1, 'extensive');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rango = calculateRangoSeleccion(P,f,presion,params) %#ok<INUSD,INUSL>
sf=length(f);
rango = zeros(sf,2);
M = max(f);
m = min(f);

if (m~=M)
    rango(1:sf,2)           = 1+presion*(((f(1:sf)'-m)/(M-m))-1);
    rango(find(isnan(f)),2) = 0; %#ok<FNDSB>
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inds = seleccion(rango, ntimes, params)
%this function selects an amount "ntimes" of individuals with a probability
%directly proportional to its fitness
%
% rango: matrix num_individualsX2, each row of the for [0..num], expressing
% the probability of selecting each individual.

if ntimes==1
  r = rango(:,1) + ( rand(size(rango,1),1).*(rango(:,2) - rango(:,1)) ); %randnumber(rango);
else
  %generates a matrix of size "size(rango,1)Xntimes" such that the i-th row
  %is an array of "ntimes" random numbers in the range [rango(i,1)...rango(i,1)]  
  r = bsxfun(@plus, rango(:,1), bsxfun(@times, rand(size(rango,1),ntimes), (rango(:,2) - rango(:,1))) );
end

%this is for the case of extremely high rates of 0 or NaN fitnesses, so the
%population is not immediately dominated by one individual. For small
%enough values of params.jiggleSelection (for example, eps(1)), this will
%become unnoticeable for high fitness
if params.jiggleSelection>0
  r = r + (rand(size(r))-0.5)*params.jiggleSelection;
end

[nevermind,inds] = max(r); %#ok<ASGLU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see the tech report "Fitness Uniform Selection to Preserve Genetic
% Diversity" from Marcus Hutter
function idxInP = FUSSfun(P, fitness, params, numsamples, samples) %#ok<INUSL>

FUSS = params.selection.FUSS;

mxf = max(fitness);
mnf = min(fitness);

%samples from the range of fitness in the population
if isempty(samples)
  samples = rand(numsamples, 1)*(mxf-mnf)+mnf;
end

dists       = abs(repmat(fitness, 1, numel(samples)) - repmat(samples', numel(fitness), 1));

idxInP  = nan(size(samples));

switch FUSS.epsMode
  case 'abs'
    epsd    = FUSS.epsDist;
  case 'rel'
    epsd    = FUSS.epsDist*(mxf-mnf);
  otherwise
    error('params.selection.FUSS.epsMode==%s not understood!!!!\n', FUSS.epsMode);
end

switch FUSS.mode
  case 'plain'
    mndists     = min(dists);
    offsetDists = dists-repmat(mndists, numel(fitness), 1)<epsd;
    for k=1:numel(samples)
      adyacent = find(offsetDists(:,k));
      if numel(adyacent)>1
        idxInP(k) = pickOption(adyacent, []);
      else
        idxInP(k) = adyacent;
      end
    end
  case 'levelled1'
    samplesAssigned   = false(size(samples)); %samples assigned to a fitness value
    enfranchisedtimes = zeros(size(fitness)); %number of times a fitness has been assigned to samples
    maxfranchise      = FUSS.ntimes;          %maximum for that number
    enfranchisedtimes(isnan(dists(:,1))) = maxfranchise; %NAN fitnesses are not elligible
    
    %get fitnesses too close between them
    tooClose          = abs(repmat(fitness, 1, numel(fitness)) - repmat(fitness', numel(fitness), 1))<epsd;
    
    distsRows         = find(~isnan(fitness)); %indexes of rows (fitnesses)
    distsCols         = 1:numel(samples); %indexes of cols (samples)
    dists             = dists(distsRows, :); %distances between rows (fitnesses) and cols (samples)
    while (~isempty(dists))
      [mnd mnd] = min(dists(:)); %get min distance
      row = mod(mnd-1, size(dists,1))+1; %fitness
      col = (mnd-row)/size(dists,1)+1;   %sample
      drow = distsRows(row); %fitness index
      dcol = distsCols(col); %sample index
      closes = find(tooClose(drow, :)); %rows (fitnesses) that are close to the minimum (that is, they are considered undistinguishable)
      if numel(closes)>1 %there are some other fitnesses close to this one. 
        drow = pickOption(closes, []); %so, select one of this group
        row  = find(distsRows==drow);  %and get its row in the dwindling matrix "dists"
        if numel(row)~=1
          error('This is very strange!!!!');
        end
        enfranchisedtimes(drow) = maxfranchise; %make sure that this fitness cannot be selected again
        %quit the selected fitness from the elligible fitnesses
        tooClose(drow,:) = false;
        tooClose(:,drow) = false;
      end
      idxInP(dcol) = drow;  %assign this fitness to this sample
      samplesAssigned(dcol) = true; %mark this sample as assigned
      dists(:,col)   = []; %get the sample out of the matrix
      distsCols(col) = [];
      enfranchisedtimes(drow) = enfranchisedtimes(drow) + 1; %add 1 to the number of times this fitness has been used
      if enfranchisedtimes(drow)>=maxfranchise %if it has been used too many times:
        dists(row,:) = []; %get the fitness out of the matrix
        distsRows(row) = [];
      end
    end
  otherwise
    error('params.selection.FUSS.mode==%s not understood!!!!\n', FUSS.mode);
end

if any(isnan(idxInP))
  isn = isnan(idxInP);
  nnew  = sum(isn);
  mprintf(params.ftst, 'evolucionGeneric:FUSS: THERE REMAIN %d SAMPLES WITHOUT ASSIGNMENT.\nTHEY WILL BE ASSIGNED BY REISUEEING A SAMPLE SELECTION\n', nnew);
  idxInP(isn) = FUSSfun(P, fitness, params, nnew, samples(isn));
end