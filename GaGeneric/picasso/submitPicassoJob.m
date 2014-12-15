function newjob = submitPicassoJob(params, terclus, newjob, numTasks, generation)

nl = sprintf('\n');

inputTemplate  = strrep(terclus.inputFileNameTemplate, '%', '%%');
outputTemplate = strrep(terclus.outputFileNameTemplate, '%', '%%');

evparams = params.evparams; %#ok<NASGU>
save([params.nomdir terclus.evparamsFileNameTemplate], 'evparams', '-v7');
clear evparams;

if terclus.showOutput
  showoutput = '> enMedio_${PBS_ARRAY_INDEX}.txt';
else
  showoutput = '';
end

if terclus.remoteSSH
  ndir = terclus.ssh.nomdir;
else
  ndir = params.nomdir;
end
ndirdeployed = ndir;

if terclus.ssh.useMatlab
  main = 'matlab -nodisplay -nojvm -r executePicassoTask"(''';
  sep = ''', ''';
  ending = [''', ''${PBS_ARRAY_INDEX}'', ''' doAddPath(terclus) ''')" '];
else
  main = ['cd deployed' nl './my_run_executePicassoTask.sh /export/home_users/home/soft/matlab/R2007b/toolbox/compiler/mcr "'];
  sep = ' ';
  ending = ' ${PBS_ARRAY_INDEX} " ';
  ndirdeployed = ['../' ndirdeployed];
end

str = [
'#!/bin/bash' nl ...
'# GENERACION ' mat2str(generation) nl ...
'# numero de cpus que empleara el calculo:' nl ...
'#PBS -l ncpus=1' nl ...
'# memoria que empleara el calculo:' nl ...
'#PBS -l mem=' terclus.maxPicassoMem nl ...
'# como mucho tardara 5 horas' nl ...
'#PBS -l walltime=' terclus.maxPicassoTime nl ...
'# si se quiere usar un array job, ahora mismo esta comentado:' nl ...
'#PBS -J 1-' mat2str(numTasks) nl ...
'# para que lo envie al superdome:' nl ...
'#PBS -q routex86' nl ...
'' nl ...
'# para que vaya al directorio actual:' nl ...
'cd "$PBS_O_WORKDIR"' nl ...
...%'cd GaGeneric/picasso' nl ...
...%'echo pasoPorAqui > pasoPorAqui_${PBS_ARRAY_INDEX}.txt' nl ...
'' nl ...
'# programa a ejecutar, con sus argumentos:' nl ...
main ndirdeployed sep inputTemplate sep outputTemplate  ending showoutput nl ...
...%'echo pasoPorAlli > pasoPorAlli_${PBS_ARRAY_INDEX}.txt' nl ...
'' nl ...
];

fname       = [params.nomdir terclus.scriptFileName];
fnameremote = [ndir terclus.scriptFileName];
f = fopen(fname, 'w');
fprintf(f, str);
fclose(f);

[FAILED, cmdOut, jobID] = sendJob(newjob, params, terclus, fname, fnameremote);
  
if FAILED
  error('qsub FAILED, got <<%s>>!!!\n', cmdOut);
end

if isempty(jobID)
  error('jobId is empty, got <<%s>> instead!!!\n', cmdOut);
end

newjob.picassoJob = true;

newjob.jobID = jobID;

newjob.numTasks = numTasks;

newjob.nomdir = params.nomdir;

newjob.inputFileNameTemplate = terclus.inputFileNameTemplate;
newjob.outputFileNameTemplate = terclus.outputFileNameTemplate;
newjob.jobIDFileName = terclus.jobIDFileName;

f = fopen([params.nomdir terclus.jobIDFileName], 'w');
fprintf(f, 'FOR NOMDIR %s, GENERATION %d, NUMTASKS %d, JOBID: %s\n', params.nomdir, generation, numTasks, mat2str(jobID));
fclose(f);

function pth = doAddPath(terclus)

jobArgs = terclus.jobArgs;

for k=1:numel(jobArgs)
  if strcmpi(jobArgs{k}, 'PathDependencies')
    pth  = jobArgs{k+1};
    pth  = [pth, repmat({pathsep}, size(pth))]; %#ok<AGROW>
    pth  = reshape(pth', 1, []);
    pth  = horzcat(pth{:});
    return
  end
end

error('This should never happen!!!!');

function [FAILED, cmdOut, jobID] = sendJob(newjob, params, terclus, fname, fnameremote)

command = ['qsub ' fnameremote];

if terclus.remoteSSH
 
  channel = createSSHChannel(params, terclus);
  
  sp = channel.createSCPClient;
  qw=javaArray('java.lang.String', numel(newjob.inputFiles)+2);
  for k=1:numel(newjob.inputFiles)
    qw(k) = java.lang.String(newjob.inputFiles{k});
  end
  qw(numel(newjob.inputFiles)+1) = java.lang.String([terclus.ssh.nomdir terclus.evparamsFileNameTemplate]);
  qw(numel(newjob.inputFiles)+2) = java.lang.String(fname);
  sp.put(qw, terclus.ssh.nomdir);
  
  [channel, cmdOut]  =  executeCommand(channel,command);
  FAILED = isempty(cmdOut) || all(isspace(cmdOut));
  
  channel.close;
  
else

  [FAILED, cmdOut] = system(command);

end

%fprintf('Submitted job with command %s, STATUS: %d, output: <<%s>>\n', command, FAILED, cmdOut);

jobID = regexp( cmdOut, '[0-9\[\]]+\.[a-zA-Z0-9-\.]*', 'match', 'once' );
