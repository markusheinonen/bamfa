%
% Plots Bayesian flux sample pairs
%
% Parameters:
% - model        : COBRA model
% - sol          : sample struct from bmfa/bfba
% - fluxpairinds : (n,2) matrix of flux pair indices
%                  e.g. [3 10; 15 16; 15 18] will plot flux pairs
%                      3-10, 15-16 and 15-18
% - fba          : fba solution to be plotted (from optimizeCbModel)
%
% Copyright (c) 2016 Markus Heinonen
%
function [] = plotfluxpair2(model, sols, fluxpairinds, fba)

	if ~iscell(sols)
		sols = {sols};
	end
	if ~exist('fluxpairinds','var')
		fluxpairinds = [1 2; 1 3; 1 4; 2 3; 2 4; 2 5];
	end
	if ~exist('fba','var')
		fba = [];
	end
	
	% 6 color "set1" color palette from colorbrewer2.org
	% red - green - blue - purple - orange - yellow
	colors = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51] ./ 255;
	alpha = 0.5;	
	
	D = size(fluxpairinds,1);
	Nj = ceil(sqrt(D));
	Ni = ceil(sqrt(D));

	for d=1:size(fluxpairinds,1)
		subplot(Ni,Nj,d);
%		i = mod(d,Ni);
%		j = mod(d,Nj);
		vi = fluxpairinds(d,1);
		vj = fluxpairinds(d,2);
		
		minx = -inf;
		maxx = inf;
		miny = -inf;
		maxy = inf;
		for k=1:length(sols)
			xi = sols{k}.vsamples(vi,:);
			xj = sols{k}.vsamples(vj,:);
			
			minx = min(minx, min(xi));
			maxx = max(maxx, max(xi));
			miny = min(miny, min(xj));
			maxy = max(maxy, max(xj));
			
	%		subplot('position',[0.08+0.90*(i-1)/Ni 0.08+0.90*(Ni-j)/Ni 0.90/Ni-0.006 0.90/Ni-0.008]);

			scatter(xi,xj, 4, colors(k,:), 'filled', 'MarkerFaceAlpha',alpha);
			hold on;
			
			if ~isempty(fba)
				plot(fba.x(vi), fba.x(vj), 'k.', 'MarkerSize',18); 
			end
		end
		hold off;

		grid;
		
		add = 0.15;
		
		xrng = maxx-minx;
		xlim([ minx-xrng*add maxx+xrng*add]);

		yrng = maxy-miny;
		ylim([ miny-yrng*add maxy+yrng*add]);

		ylabel(model.rxns(vj), 'interpreter','none');
		xlabel(model.rxns(vi), 'interpreter','none');
	end
end


