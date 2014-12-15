function doEvoRuns(modo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
args_cluster_pica1 = {
    'tam_poblacion', '100', ...
    'generaciones', '100', ...
    'terclus.concurrentWorkers', '100', ...
    'terclus.pauseTime', '60', ...
    'terclus.pauseTimeAfterError', '300', ...
    'terclus.maxPicassoTime', '''1:00:00''', ...
    'terclus.remoteSSH', 'true', ...
    'terclus.ssh.password', '''f2367c21''', ...
  };
args_cluster_pica2 = {
    'tam_poblacion', '100', ...
    'generaciones', '150', ...
    'terclus.concurrentWorkers', '4', ...
    'terclus.pauseTime', '600', ...
    'terclus.pauseTimeAfterError', '600', ...
    'terclus.maxPicassoTime', '''5:00:00''', ...
    'terclus.remoteSSH', 'true', ...
    'terclus.ssh.password', '''f2367c21''', ...
  };
args_cluster_pica21 = {
    'tam_poblacion', '100', ...
    'generaciones', '150', ...
    'terclus.concurrentWorkers', '4', ...
    'terclus.pauseTime', '1200', ...
    'terclus.pauseTimeAfterError', '1200', ...
    'terclus.maxPicassoTime', '''5:00:00''', ...
    'terclus.remoteSSH', 'true', ...
    'terclus.ssh.password', '''f2367c21''', ...
  };
args_cluster_pica3 = {
    'tam_poblacion', '100', ...
    'generaciones', '200', ...
    'terclus.concurrentWorkers', '1', ...
    'terclus.pauseTime', '1800', ...
    'terclus.pauseTimeAfterError', '1800', ...
    'terclus.maxPicassoTime', '''9:00:00''', ...
    'terclus.remoteSSH', 'true', ...
    'terclus.ssh.password', '''f2367c21''', ...
  };
args_cluster_pica4 = {
    'tam_poblacion', '100', ...
    'generaciones', '100', ...
    'terclus.concurrentWorkers', '4', ...
    'terclus.pauseTime', '600', ...
    'terclus.pauseTimeAfterError', '1800', ...
    'terclus.maxPicassoTime', '''5:00:00''', ...
    'terclus.remoteSSH', 'true', ...
    'terclus.ssh.password', '''f2367c21''', ...
  };
args_cluster_pica5 = {
    'tam_poblacion', '100', ...
    'generaciones', '200', ...
    'terclus.concurrentWorkers', '4', ...
    'terclus.pauseTime', '600', ...
    'terclus.pauseTimeAfterError', '1800', ...
    'terclus.maxPicassoTime', '''5:00:00''', ...
    'terclus.remoteSSH', 'true', ...
    'terclus.ssh.password', '''f2367c21''', ...
  };
args_cluster_terc1 = {
    'tam_poblacion', '100', ...
    'generaciones', '100', ...
    'terclus.concurrentWorkers', '100', ...
    'terclus.pauseTime', '1', ...
    'terclus.pauseTimeAfterError', '60', ...
  };
args_cluster_terc2 = {
    'tam_poblacion', '100', ...
    'generaciones', '150', ...
    'terclus.concurrentWorkers', '100', ...
    'terclus.pauseTime', '1', ...
    'terclus.pauseTimeAfterError', '60', ...
  };
args_cluster_terc3 = {
    'tam_poblacion', '500', ...
    'generaciones', '500', ...
    'terclus.concurrentWorkers', '100', ...
    'terclus.pauseTime', '1', ...
    'terclus.pauseTimeAfterError', '60', ...
  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
args_mer_inchworm = {
    'evparams.fitnessCalc', '''mode3''', ...
    'evparams.walker.row.attachMode', '1', ...
    'evparams.walker.d3.borderRating', '''method1''', ...
    'evparams.walker.d3.rotationMode', '''inchworm''' ...
  };

args_mer_handoverhand = {
    'evparams.fitnessCalc', '''hoh''', ...
    'evparams.walker.row.attachMode', '0', ...
    'evparams.walker.d3.borderRating', '''handoverhand1''', ...
    'evparams.walker.d3.rotationMode', '''handoverhand''' ...
  };

args_mer_handoverhand2 = {
    'evparams.fitnessCalc', '''mode3''', ...
    'evparams.walker.row.attachMode', '0', ...
    'evparams.walker.d3.borderRating', '''handoverhand1''', ...
    'evparams.walker.d3.rotationMode', '''handoverhand''' ...
  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
args_sel_paco = {
      'selection.mode', '''PACO''', ...
      'evparams.fitness.HOH.zeroIsNan', 'false', ...
      'rewindNans', 'false', ...
      'nansReplaced', 'false', ...
      'elitism', 'true' ...
  };

args_sel_fuss = {
      'selection.mode', '''FUSS''', ...
      'evparams.fitness.HOH.zeroIsNan', 'true', ...
      'rewindNans', 'true', ...
      'nansReplaced', 'true', ...
      'elitism', 'false' ...
  };

switch modo
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1
    conf         = molMot3DConf(false);
    args_cluster = args_cluster_pica2;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_inchworm;
    subdir       = ['resultados' filesep 'INCHWORM_PACO'];
  case 2
    conf         = molMot3DConf(true);
    args_cluster = args_cluster_pica2;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_inchworm;
    subdir       = ['resultados' filesep 'INCHWORM_FO_PACO'];
  case 3
    conf         = molMot3DConf(false);
    args_cluster = args_cluster_pica21;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_handoverhand2;
    subdir       = ['resultados' filesep 'HANDOVERHAND_PACO'];
  case 4
    conf         = molMot3DConf(true);
    args_cluster = args_cluster_pica21;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_handoverhand2;
    subdir       = ['resultados' filesep 'HANDOVERHAND_FO_PACO'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 5
    conf         = molMot3DConf(false);
    args_cluster = args_cluster_terc2;
    mode         = 'remote';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_inchworm;
    subdir       = ['resultados' filesep 'INCHWORM_PACO'];
  case 6
    conf         = molMot3DConf(true);
    args_cluster = args_cluster_terc2;
    mode         = 'remote';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_inchworm;
    subdir       = ['resultados' filesep 'INCHWORM_FO_PACO'];
  case 7
    conf         = molMot3DConf(false);
    args_cluster = args_cluster_terc2;
    mode         = 'remote';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_handoverhand2;
    subdir       = ['resultados' filesep 'HANDOVERHAND_PACO'];
  case 8
    conf         = molMot3DConf(true);
    args_cluster = args_cluster_terc2;
    mode         = 'remote';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_handoverhand2;
    subdir       = ['resultados' filesep 'HANDOVERHAND_FO_PACO'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 9
    conf         = molMot3DConf(false);
    args_cluster = args_cluster_pica3;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_inchworm;
    subdir       = ['resultados' filesep 'INCHWORM_PACO'];
  case 10
    conf         = molMot3DConf(true);
    args_cluster = args_cluster_pica3;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_inchworm;
    subdir       = ['resultados' filesep 'INCHWORM_FO_PACO'];
  case 11
    conf         = molMot3DConf(false);
    args_cluster = args_cluster_pica3;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_handoverhand2;
    subdir       = ['resultados' filesep 'HANDOVERHAND_PACO'];
  case 12
    conf         = molMot3DConf(true);
    args_cluster = args_cluster_pica3;
    mode         = 'picasso';
    args_sel     = args_sel_paco;
    args_mer     = args_mer_handoverhand2;
    subdir       = ['resultados' filesep 'HANDOVERHAND_FO_PACO'];
end


for k=1:10;
  conf.samples.mainfun(mode, {subdir 'newsim'}, ...
    'evparams.fitness.HOH.mode', '''absabs''', ...
    'newsQuotaMode', '''custom''', ...
    'evparams.walker.d3.hoh.swingAngle', 'pi/18', ...
    'newsQuota', '@(gen, rid, genInRange, range, P, params) ((sum(isnan(P.fitness))==numel(P.fitness)) * numel(P.fitness)) + ((sum(isnan(P.fitness))~=numel(P.fitness)) * round(sum(isnan(P.fitness))*0.75))', ...
    'genome.poisson', '{0.5, 0:10, poisspdf(0:10, 0.5)}', ...
    'evparams.walker.row.spchg.constT', '2', ...
    args_cluster{:}, ...
    args_sel{:}, ...
    args_mer{:} ...
  );
end
