%
% Plots Bayesian flux sample pairs on grid
%
% Parameters:
% - model    : COBRA model
% - sol      : sample struct from bmfa/bfba
% - fluxinds : flux indices to be plotted
%              warning: more than 10 fluxes do not draw well in matlab
% - fba      : fba solution to be plotted (from optimizeCbModel)
%
% Copyright (c) 2016 Markus Heinonen
%
function [] = plotfluxpair(model, sol, fluxinds, fba)

	if ~exist('fluxinds','var')
		target = find(model.c);
		target = target(1:min(3,length(target)));
		[~,top] = sort(std( sol.vsamples(:,:)'), 'descend');
		fluxinds = [target top(1:5)];
	end
	if ~exist('fba','var')
		fba = [];
	end
	
	% 6 color "set1" color palette from colorbrewer2.org
	% red - green - blue - purple - orange - yellow
	colors = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 255 255 51] ./ 255;
	alpha = 0.5;

	I = fluxinds;
	Np = 50;    % points per density axis
	Ni = length(I);
	Ns = size(sol.vsamples(:,:),2);
	
	if Ni > 20
		warning('Too many fluxes, maximum number of fluxes is 20');
		return;
	end
	
	[N,K,C] = size(sol.vsamples);
	Vs = reshape(sol.vsamples, N, K*C);
		
	for i=1:Ni
		for j=i:Ni
			
%			subplot(Ni,Ni,Ni*(j-1)+i);
			subplot('position',[0.08+0.90*(i-1)/Ni 0.08+0.90*(Ni-j)/Ni 0.90/Ni-0.006 0.90/Ni-0.008]);
			%                   left               bottom              width   height
			
			vi = I(i);
			vj = I(j);

			minvi = min(Vs(vi,:));
			maxvi = max(Vs(vi,:));
			minvj = min(Vs(vj,:));
			maxvj = max(Vs(vj,:));
			
			minvi = min(minvi, minvi-(0.03-(maxvi-minvi))/2);
			maxvi = max(maxvi, maxvi+(0.03-(maxvi-minvi))/2);
			minvj = min(minvj, minvj-(0.03-(maxvj-minvj))/2);
			maxvj = max(maxvj, maxvj+(0.03-(maxvj-minvj))/2);
			
			% add 6% to both sides
			rngi = (maxvi-minvi);
			rngj = (maxvj-minvj);
			minvi = minvi - 0.05*rngi;
			maxvi = maxvi + 0.05*rngi;
			minvj = minvj - 0.05*rngj;
			maxvj = maxvj + 0.05*rngj;
				
			if vi~=vj
%				X = Vs([vi vj],:)';
				
				xi = linspace(min(Vs(vi,:))-0.01, max(Vs(vi,:))+0.01, Np)';
				xj = linspace(min(Vs(vj,:))-0.01, max(Vs(vj,:))+0.01, Np)';
				[x11,x22] = meshgrid(xi,xj);
				XI = [x11(:) x22(:)];
				F = reshape(ksdensity(Vs([vi vj],:)',XI),Np,Np);
		
				maxval = max(F(:));
				
				if Ns < 1000
					[c,h] = contourf(xi,xj,F, [0 linspace(maxval/20,maxval,4)], 'edgecolor','none');
					hold on;
				end
				scatter(Vs(vi,:), Vs(vj,:), 4, colors(1,:), 'filled', 'MarkerFaceAlpha',alpha);
				hold on;
%				plot(sol.vmap(vi), sol.vmap(vj), 'r.', 'MarkerSize',18);
				
				if ~isempty(fba)
					plot(fba.x(vi), fba.x(vj), 'k.', 'MarkerSize',18); 
				end
				
				if model.lb(vi) > -1000
					fill([minvi model.lb(vi) model.lb(vi) minvi], [minvj minvj maxvj maxvj], 'black','facealpha',0.2,'edgecolor','none');
				end
				if model.lb(vj) > -1000
					fill([minvi maxvi maxvi minvi], [minvj minvj model.lb(vj) model.lb(vj)], 'black','facealpha',0.2,'edgecolor','none');
				end
				
				hold off;

				caxis([-maxval maxval]);
				colormap(redblue);
				xlim([minvi maxvi]);
				ylim([minvj maxvj]);
			else
				histogram(Vs(vi,:))
				xlim([minvi maxvi]);
			end
			
			if i~=j
				grid;
			end
			
			if i==1 && j==Ni
				ylabel(model.rxns(vj), 'interpreter','none');
				xlabel(model.rxns(vi), 'interpreter','none');
			elseif i==1
				ylabel(model.rxns(vj), 'interpreter','none');
				set(gca,'xticklabel',[]);
			elseif j==Ni
				xlabel(model.rxns(vi), 'interpreter','none');
				set(gca,'yticklabel',[]);
			else
				set(gca,'yticklabel',[],'xticklabel',[]);
			end
			
			if i==1 && j==1
				set(gca,'yticklabel',[]);
			end
		end
	end
end


