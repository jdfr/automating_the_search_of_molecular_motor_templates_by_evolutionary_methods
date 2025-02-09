function experiments = reapAllStats(basedir)
%script to reap all statistics

experiments = { ...
  ... %EN MI ORDENADOR
01, 'mios\03_Apr_2009_20_52_03\1_0.0001_1', 1000, 1/3; ... % no muy v�lido (n�mero de Gs anormalmente bajo)
02, 'mios\05_Apr_2009_19_47_11\1=0.75_0.0001_1', 76, 0.75; ...% max G=50
03, 'mios\05_Apr_2009_21_29_05\1=0.75_0.0001_1', 141, 0.75; ...
04, 'mios\06_Apr_2009_00_34_53\1=0.75_0.0001_1', 1000, 0.75; ...
05, 'mios\06_Apr_2009_00_34_53\1=0.75_0.0001_2', 1000, 0.75; ...
06, 'mios\06_Apr_2009_00_34_53\1=0.75_0.0001_3', 1000, 0.75; ...
07, 'mios\08_Apr_2009_18_27_48\1=0.5_0.0001_1', 1000, 0.5; ...
08, 'mios\08_Apr_2009_18_27_48\1=0.5_0.0001_2', 1000, 0.5; ...
09, 'mios\08_Apr_2009_18_27_48\1=0.5_0.0001_3', 250, 0.5; ...
10, 'mios\09_Apr_2009_11_53_49\1=0.33333_0.0001_1', 1000, 1/3; ...
11, 'mios\09_Apr_2009_11_53_49\1=0.33333_0.0001_2', 1000, 1/3; ...
12, 'mios\09_Apr_2009_11_53_49\1=0.33333_0.0001_3', 1000, 1/3; ...
13, 'mios\09_Apr_2009_11_53_49\1=0.33333_0.0001_4', 1000, 1/3; ...
14, 'mios\09_Apr_2009_11_53_49\1=0.33333_0.0001_5', 1000, 1/3; ...
15, 'mios\09_Apr_2009_11_53_49\2=0.2_0.0001_1', 1000, 0.2; ...
16, 'mios\09_Apr_2009_11_53_49\2=0.2_0.0001_2', 1000, 0.2; ...
17, 'mios\09_Apr_2009_11_53_49\2=0.2_0.0001_3', 196, 0.2; ...
18, 'mios\11_Apr_2009_13_18_00\1=0.2_0.0001_1', 356, 0.2; ...
19, 'mios\11_Apr_2009_17_04_47\1=0.2_0.0001_1', 1000, 0.2; ...
20, 'mios\11_Apr_2009_17_04_47\1=0.2_0.0001_2', 690, 0.2; ...
21, 'mios\12_Apr_2009_10_44_51\1=1_0.0001_1', 1000, 1; ...
22, 'mios\12_Apr_2009_10_44_51\1=1_0.0001_2', 1000, 1; ...
23, 'mios\12_Apr_2009_10_44_51\1=1_0.0001_3', 1000, 1; ...
24, 'mios\12_Apr_2009_10_44_51\1=1_0.0001_4', 413, 1; ...
  ... %EN EL DE GEMA
25, 'gema\09_Apr_2009_12_27_50\1=0.0001_0.0001_1', 1000, 0.0001; ...
26, 'gema\09_Apr_2009_12_27_50\1=0.0001_0.0001_2', 1000, 0.0001; ...
27, 'gema\09_Apr_2009_12_27_50\1=0.0001_0.0001_3', 1000, 0.0001; ...
28, 'gema\09_Apr_2009_12_27_50\1=0.0001_0.0001_4', 1000, 0.0001; ...
29, 'gema\09_Apr_2009_12_27_50\1=0.0001_0.0001_5', 1000, 0.0001; ...
30, 'gema\09_Apr_2009_12_27_50\2=0.1_0.0001_1', 613, 0.1; ...
31, 'gema\12_Apr_2009_10_55_43\1=0.1_0.0001_1', 1000, 0.1; ...
32, 'gema\12_Apr_2009_10_55_43\1=0.1_0.0001_2', 807, 0.1 ...
};

nranges = 2;

experiments = struct('path', experiments(:,2), 'generations', experiments(:,3), 'presion', experiments(:,4), 'info', cell(size(experiments,1),1));

for k=1:numel(experiments)
  fprintf('Reaping stats from %d generations for %s...\n', experiments(k).generations, experiments(k).path);
  experiments(k).info = reapTimeStats([basedir filesep experiments(k).path],experiments(k).generations,nranges,false);
end
fprintf('all done!\n');
