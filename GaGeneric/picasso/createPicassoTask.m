function newjob = createPicassoTask(params, terclus, numTask, generation, newjob, evfun, noutputs, inputs)
task.Function                = evfun;
task.NumberOfOutputArguments = noutputs;
task.InputArguments          = inputs;
task.numTask                 = numTask;
task.generation              = generation;
task.inputFile               = sprintf(terclus.inputFileNameTemplate, numTask);
if terclus.remoteSSH
  ndir = terclus.ssh.nomdir;
else
  ndir = params.nomdir;
end
newjob.inputFiles{numTask}   = [ndir task.inputFile];
newjob.outputFiles{numTask}  = [ndir sprintf(terclus.outputFileNameTemplate, numTask)];
save([params.nomdir task.inputFile], 'task', '-v7');

