%
% Plots Bayesian flux covariances
%
% Parameters:
% - model    : COBRA model
% - sols     : sample struct from bmfa/bfba (a cell of multiple sols can be given)
% - fluxinds : flux indices
%
% Copyright (c) 2016 Markus Heinonen
%
function [] = plotfluxcov(model, sol, fluxinds)
	
	N = size(model.S,2);
	
	if ~exist('fluxinds','var')
		fluxinds = 1:min(N,50);
	end
	Nids = nnz(fluxinds);


	[N,Ns,Nc] = size(sol.vsamples);

	X = sol.vsamples(fluxinds,:,:);
	X = reshape(X, [Nids Ns*Nc])';
	
	C = corr(X);

	% order columns by first graph laplacian component
	L = diag(sum(C)) - C;
	[coeff,score,latent] = pca(L);
	[v,p] = sort(score(:,1));

	% go through all columns, swap closest column to be next
	for i=1:Nids-1
%		[~,k] = min(sum(( A(:,i) - A(:,i+1:end)).^2));		
		[~,k] = min(sum((  C(p,p(i)) - C(p, p(i+1:end))).^2));		
		
		% swap i+1 and k+i
		t = p(i+1);
		p(i+1) = p(i+k);
		p(i+k) = t;
	end
	
	Cpp = C(p,p);
%	Cpp = tril(Cpp);
	
	imagesc(Cpp);
	
	hold on;
	for i = 1:95
	   plot([.5,76.5],[i-.5,i-.5],'k-', 'color', 0.8*[1 1 1]);
	   plot([i-.5,i-.5],[.5,76.5],'k-', 'color', 0.8*[1 1 1]);
	end
	hold off;
	
	caxis([-1 1]);
	colorbar;
	colormap(redbluecmap);

	title('Flux correlations');
	
	set(gca,'TickLabelInterpreter','none');
	set(gca,'yticklabel',model.rxns(fluxinds(p)),'ytick',1:Nids, 'Ticklength', [0 0]);
	set(gca,'xticklabel',model.rxns(fluxinds(p)),'xtick',1:Nids, 'XTickLabelRotation',90);
	
	
	% exclude [-0.1 0.1]
	A = C(p,p);
	A(abs(A)<0.1) = 0;

end

