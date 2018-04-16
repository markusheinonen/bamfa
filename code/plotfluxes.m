%
% Plots Bayesian flux sample distributions
%
% Parameters:
% - model    : COBRA model
% - sols     : sample struct from bmfa/bfba (a cell of multiple sols can be given)
% - fluxinds : flux indices
% - fba      : fba solution to be plotted (from optimizeCbModel)
%
% Copyright (c) 2016 Markus Heinonen
%
function [] = plotfluxes(model, sols, fluxinds, fba)
	
	if ~iscell(sols)
		sols = {sols};
	end
	if ~exist('fluxinds', 'var')
		fluxinds = 1:10;
	end
	if ~exist('fba', 'var')
		fba = [];
	end

	N = nnz(fluxinds);
	rows = ceil(sqrt(N));
	cols = ceil(sqrt(N));
	
	% 6 color "set1" color palette from colorbrewer2.org
	% red - green - blue - purple - orange - yellow
	colors = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51] ./ 255;
	alpha = 0.5;
	
	lb = model.lb;
	ub = model.ub;
	
	for i=1:length(fluxinds)
		subplot(rows,cols*2,[i*2-1 i*2]);

%		r = rem(i-1,6)+1;
%		r = floor(i/13)+1;
%		s = rem(i-1,6)+1;
%		subplot('position',[0.08+0.90*(s-1)/Nj 0.08+0.90*(Ni-r)/Ni 0.90/Nj-0.004 0.90/Ni-0.055]);

		j = fluxinds(i);
		
		hold on;
		minx = inf; maxx = -inf;
		for k=1:length(sols)
			y = sols{k}.vsamples(j,:);

			[f,x] = ksdensity(y);
			f = f / max(f);
			Ibnd = (x >= lb(j)) & (x <= ub(j)); % clip at bounds
			f(~Ibnd) = 0;

			I = f > 0.01; % clip tiny values

			minx = min(minx, min(x(I)));
			maxx = max(maxx, max(x(I)));

			area(x(I),f(I), 'facecolor', colors(k,:), 'facealpha',alpha, 'edgecolor', colors(k,:), 'edgealpha',alpha, 'linewidth',1.5);
		end
		
		if ~isempty(fba)
			plot(fba.x(j), 0, 'k.', 'MarkerSize',18); 
		end
		
		hold off;
%		rangex = maxx-minx;
%		if nnz(Ibnd)
%			xlim([minx-0.04*rangex-0.01 maxx+0.04*rangex+0.01]);
%			hold on;
%			plot([minx-0.04*rangex-0.01 maxx+0.04*rangex+0.01], [0 0], 'k-', 'linewidth',1.5);
%			hold off;
%		end

%		set(gca, 'Ytick', [], 'YTickLabel', [], 'box','off','YColor','w');
		title(model.rxns(j), 'Interpreter','none');
		set(gca, 'color','w');
		xlim('auto');
		ylim([0 1]);
		grid on;
	end
end




