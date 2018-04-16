function [sol,model] = bmfa(model, y, g, sdx, Nsamples, Nchains, sampler, T, svdtol)
% Bayesian Metabolic Flux Analysis
% 
% Given a set of known flux values (e.g. inputs/outputs) infers the
%   distribution of the unknown fluxes. Upholds steady-state with error
%   'sdx' (standard deviation)
%
% Inputs:
% - model    : COBRA model, a struct with fields {S,lb,ub}
% - y        : a (N x 1) sparse vector of known fluxes
% - g        : a (N x 1) sparse vector of known flux errors (std)
% - sdx      : steady-state error (std) [constant or (N x 1) vector]
% - Nsamples : number of samples to simulate [default 200]
% - Nchains  : number of chains to simulate  [default 5]
% - sampler  : 'gibbs' [default] or 'hmc' 
% - T        : sampling parameter, skip every T'th sample for gibbs [default 20],
%               running time per sample for hmc [default pi/2]
% - svdtol   : SVD simplication of the stoichiometry threshold [default 0]
%
% Outputs:
% struct with fields
%  .vsamples  : flux samples [ Nfluxes x Nsamples x Nchains ]
%  .dxsamples : sample dx's
%  .logps     : sample log densities (higher indicates better solution)
%  .vmu       : flux posterior mean
%  .vcov      : flux posterior cov
%  .vmap      : best flux solution
%  .r         : potential scale reduction factors of fluxes, below 1.1 indicates
%               sampling convergences. Higher values mean that some rare
%               flux combinations were not found
%
% Notes:
%  - 'sdx'    determines how close the results are to steady-state. For
%             numerical stability do not decrease to zero
%  - 'T'      denotes how well the flux space is sampled, increasing this
%             gives better estimate of the flux distributions, but
%             increases running time
%
% Copyright (c) 2017 Markus Heinonen, Aalto University, Finland
%

	% read from file if filename
	if ischar(model)
		model = readmodel(model);
	end
	
	tol = 0.0001;
	S = model.S;
	lb = model.lb;
	ub = model.ub;
	[M,N] = size(S);

	if ~exist('y', 'var') || isempty(y)
		y = sparse(N,1);
	end
	if ~exist('g', 'var')
		g = sparse(N,1);
	end
	if ~exist('sdx', 'var')
		sdx = 0.01;
	end
	if ~exist('Nsamples', 'var')
		Nsamples = 100;
	end
	if ~exist('Nchains', 'var')
		Nchains = 10;
	end
	if ~exist('sampler', 'var')  % set to gibbs for over 200 reactions
		sampler = 'gibbs';
	end
	if ~exist('T', 'var')
		if strcmp(sampler,'gibbs')
			T = 50;
		else
			T = pi/2;
		end
	end
	if ~exist('svdtol', 'var')
		svdtol = 0;
	end	
	
	% change 1000's into infs so that they are ignored during sampling
	ub(ub >= 1000) = Inf;
	lb(lb <= -1000) = -Inf;

	% parameter handling
	if length(sdx) == 1 
		sdx = sdx * ones(M,1);
	end
	
	Iobs = find(g);
	G = diag(g(Iobs).^2);
	
	% means
	mdx = zeros(M,1);
	mv = min(max(0,lb),ub); % ensure bounds

	% optimise flux prior variances by MLL
	fprintf('Optimising flux variances...\n');
	sv0 = 100*ones(N,1);
	sv = sv0;
%	options = optimoptions('fmincon','display','none', 'MaxFunEvals',70000);
%	sv = fmincon(@(sv) -grad_ard(y,g,mv,sv,sdx,S),sv0,[],[],[],[],zeros(N,1),[],[],options);
%	sv = sv*2; % account for bounds

%	x0 = [mv0;sv0];
%	x = fmincon(@(x) -grad_ard(y,g,x(1:N),x(N+1:N*2),sdx,S),x0+0.1,[],[],[],[],[model.lb; zeros(N,1)],[model.ub, inf(N,1)],[],options);
%	mv = x(1:N);
%	sv = x(N+1:N*2);

	% variances
	sv(lb==ub) = tol; % zero variance if lb(i)==ub(i)
	
	% covariances
	Sv = diag(sv.^2);
	Sdx = diag(sdx.^2);
	
	% flux model
	A = (Sv*S') / (S*Sv*S' + tol*eye(M));
	vmu = mv + A*(mdx - S*mv);
	vC = Sv - A*S*Sv + A*Sdx*A';
	vC = 0.5*(vC + vC'); % ensure symmetry
	vC = vC + 1e-8*eye(N); % ensure SDP
	
	% flux posterior
	B = vC(:,Iobs) / (vC(Iobs,Iobs) + G);
	vmuy = vmu + B * (y(Iobs) - vmu(Iobs));
	vCy = vC - B * vC(Iobs,:);
	vCy = 0.5*(vCy + vCy'); % ensure symmetry
		
	% only sample constrained and non-zero variance dimensions
	Ivar = diag(vC) ~= 0;
	Ivar2 = lb~=ub;
	Ifin = ~(isinf(lb) .* isinf(ub));
	I = logical(Ivar .* Ifin .* Ivar2);	
	
	% MAP solution
	tic;
	vmap = tmvn_map(vmuy, vCy, lb, ub);
	
%	save('tmvn_case4.mat','vmuy','vCy','lb','ub','vmap');
	
	vsamples = nan(N, Nsamples, Nchains); % D x N x C
	switch sampler
		case 'hmc'
			if Nchains == 1
				vsamples(I,:) = tmvn_hmc(vmuy(I), vCy(I,I), lb(I), ub(I), vmap(I), Nsamples, T, svdtol);
			else
				vsamples(I,:,:) = tmvn_hmc2(vmuy(I), vCy(I,I), lb(I), ub(I), vmap(I), Nsamples, Nchains, T, svdtol);
			end
		case 'gibbs'
			vsamples(I,:,:) = tmvn_gibbs2(vmuy(I), vCy(I,I), lb(I), ub(I), vmap(I), Nsamples, Nchains, T, svdtol);
		case 'gpu'
			[vsamples(I,:),vmap(I)] = HMC_gpu(vmu(I), vC(I,I), lb(I), ub(I), Nsamples, T, svdtol);
	end
	runtime = toc;
	
	% extrapolate free fluxes from sampled truncated fluxes
	disp('Estimating unbounded fluxes..');
	vsamples = conditional_mvn(vsamples, vmuy, vCy, I);
%	vmap = conditional_mvn(vmap, vmuy, vCy, I);

	% compute log-densities
	disp('Computing sample log probabilities..');
	V = reshape(vsamples, [N Nsamples*Nchains]);
	ps = logmvnpdf(V, vmuy, vCy);
	logps = reshape(ps, Nsamples, Nchains);
	
	% retrieve growth and changes
	disp('Computing sample mass balances..');
	dxsamples = cell2mat( cellfun(@(X) S*X, num2cell(vsamples, [1 2]), 'un',0) );
	
	sol.sv = sv;
	sol.vmap = vmap;
	sol.vmu = vmuy;
	sol.vcov = vCy;
	sol.vmap = vmap;
	sol.vsamples = vsamples;
	sol.dxsamples = dxsamples;
	sol.logps = logps;
	sol.y = y;
	sol.g = g;
	sol.sdx = sdx;
	
	[sol.r,sol.neff] = psrf( permute(sol.vsamples, [2 1 3]) );
	sol.runtime = runtime;
	sol.Ibnd = I;
%	sol.neff = neff;
end


