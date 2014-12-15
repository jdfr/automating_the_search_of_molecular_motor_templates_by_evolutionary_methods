function showStats(experiments, indexes)

gens = [experiments(indexes).generations];

f = figure;

x=[];

for k=1:numel(indexes)
  if numel(x)~=gens(k)
    x = 1:gens(k);
  end
  subplot(numel(indexes),1,k);
  maxl = max(max([experiments(indexes(k)).info(:,1).maxl]));
  data = {...
       ...%{x, [experiments(indexes(k)).info(:,1).ming],                           'Color', 'r', 'min g'}, ...
       {x, [experiments(indexes(k)).info(:,1).meang],                          'Color', 'r', 'mean g'}, ...
       {x, [experiments(indexes(k)).info(:,1).maxg],                           'Color', 'r', 'max g'}, ...
       ...%{x, [experiments(indexes(k)).info(:,1).minl],                           'Color', 'c', 'min l'}, ...
       {x, [experiments(indexes(k)).info(:,1).meanl],                          'Color', 'c', 'mean l'}, ...
       {x, [experiments(indexes(k)).info(:,1).maxl],                           'Color', 'c', 'max l'}, ...
       ...%{x, diff([0, [experiments(indexes(k)).info(:,1).totalClusterSpawned]]), 'Color', 'm', 'time by gen.'}, ...
       {x, [experiments(indexes(k)).info(:,1).bestFitness]*maxl*2,               'Color', 'b', 'LineWidth', 2, 'best Fitness'}, ...
       };
  hold on;
  tolegend = cellfun(@(x)x{end}, data, 'uniformoutput', false);
  for m=1:numel(data)
    line(data{m}{1:end-1});
  end
  set(gca, 'XLim', [1 max(gens)]);
  grid on;
  if k==1
    legend(tolegend);
  end
  title(sprintf('presion: %.3f, comp. time: %10f', experiments(indexes(k)).presion, max(max([experiments(indexes(k)).info.totalClusterReaped])) ));
end