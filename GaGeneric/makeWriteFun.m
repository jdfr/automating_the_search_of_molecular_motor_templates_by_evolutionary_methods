function fun = makeWriteFun(params)
indivStruct = params.indivStructFull;
% conf.params.statNamesFull           = {'fitness', 'rndSeedDev', 'rndSeedEv', 'devTimeSpent', 'evTimeSpent', 'cpuTime', 'ncells', 'ncellsD', 'ncellsTeo', 'npoints', 'nsprings1', 'nsprings2', 'developed', 'overload1', 'overload2', 'atpChg1', 'atpChg2', 'offsetAbs', 'offsetRel'  ; ...
%                                        nan,       nan,          nan,         nan,            nan,           nan,       nan,      nan,       nan,         nan,       nan,         nan,         false,       nan,         nan,         nan,       nan,       nan,         nan          ; ...
%                                        '%22s',    '%s',         '%s',        '%-06.3f',      '%-06.3f',     '%-06.3f', '%03d',   '%03d',    '%03d',      '%03d',    '%03d',      '%03d',      '%d',        '%13d',      '%13d',      '%02d',    '02%d',    '%13d',      '%13d'       ;...
%                                        01,        06,           07,          03,             04,            05,        08,       09,        10,          11,        12,          13,          14,          15,          16,          17,        18,        19,          20,          ; ...
%                                        @mat2str,  @hex2num,     @hex2num,    '',             '',            '',        '',       '',        '',          '',        '',          '',          '',          '',          '',          '',        '',        '',          ''           ; ...
%                                        '%f',      '%s',         '%s',        '%f',           '%f',          '%f',      '%f',     '%f',      '%f',        '%f',      '%f',        '%f',        '%f',        '%f',        '%f',        '%f',      '%f',      '%f',        '%f'        };                                                                                      

names   = indivStruct(1,:);
chars   = indivStruct(3,:);
order   = cell2arrayMX(indivStruct(4,:));
funcs   = indivStruct(5,:);
written = cellfun('prodofsize', chars)>0;
indexes = uint8(0:(numel(names)-1))';

[order order] = sort(order);

names   =   names(order);
chars   =   chars(order);
funcs   =   funcs(order);
written = written(order);
indexes = indexes(order);

names   =   names(written);
chars   =   chars(written);
funcs   =   funcs(written);
indexes = indexes(written); %this would be required on a MEX implementation of makeParams and/or writeFun

chars = [{'%03d %d %03d'}, chars; {' '}, array2cell(repmat(' ', size(chars)))];
chars{end} = '\n';
chars = horzcat(chars{:});

fun         = writeFunFunFun(chars, names, funcs);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function fun = writeFunFunFunMX(str, indexes, funcs)
% fun = @(file, varargin)writeFunFunMX(str, indexes, funcs, file);
% 
% function fun = writeFunFunMX(str, indexes, funcs, file)
% %fun = @(gen, rid, idx, pop)writeFunMX(str, indexes, funcs, file, gen, rid, idx, pop);
% fun = @(gen, rid, idx, pop)writeFunMXCaller(str, indexes, funcs, file, gen, rid, idx, pop);
% 
% function writeFunMXCaller(str, indexes, funcs, file, gen, rid, idx, pop)
% params = writeFunMX(indexes, funcs, idx, pop);
% fprintf(file, str, gen, rid, idx, params{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fun = writeFunFunFun(str, names, funcs)
fun = @(file, varargin)writeFunFun(str, names, funcs, file);

function fun = writeFunFun(str, names, funcs, file)
fun = @(gen, rid, idx, pop)writeFun(str, names, funcs, file, gen, rid, idx, pop);
%fun = @(gen, rid, idx, pop)writeFunMX(str, names, funcs, file, gen, rid, idx, pop);

function writeFun(str, names, funcs, file, gen, rid, idx, pop)
params = makeParams(names, funcs, idx, pop);
fprintf(file, str, gen, rid, idx, params{:});

function params = makeParams(names, funcs, idx, P)
params = cell(size(names));

for k=1:numel(params)
  if isempty(funcs{k})
    params{k} = P.(names{k})(idx);
  else
    params{k} = feval(funcs{k}, (P.(names{k})(idx)));
  end
end