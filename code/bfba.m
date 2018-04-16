function [sol,model] = bfba(model, fba, sdx, Nsamples, Nchains, sampler, T, svdtol)
% Bayesian FBA
% Return a sample from the joint FBA flux distribution where:
%     v ~ U(lb,inf)
%     v ~ N(0,  sigmav.^2) 
%    Sv ~ N(0, sigmadx.^2)
%   c'v ~ N(1,  sigmac^2)
%
% Inputs:
% - model   : COBRA model, a struct with fields {S,lb,ub,c} necessary
% - fba     : fba solution [from COBRA's optimizeCbModel()]
%
% Outputs:
% struct with fields
%  .vsamples  : flux samples
%  .dxsamples : sample dx's
%  .csamples  : sample growths
%  .logps     : sample log densities (higher indicates better solution)
%  .vmu       : flux posterior mean
%  .vcov      : flux posterior cov
%
% Copyright (c) 2017 Markus Heinonen
	
	if ischar(model)
		model = readmodel(model);
	end
	if ~exist('sdx', 'var')
		sdx = 0.01;
	end
	if ~exist('Nsamples', 'var')
		Nsamples = 200;
	end
	if ~exist('Nchains', 'var')
		Nchains = 50;
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

	% fluxes
	N = size(model.S,2);
	target = find(model.c);

	% specify known fluxes as sparse vectors
	y = sparse(N,1);
	g = sparse(N,1);

	y(target) = fba.x(target);         % specify growth target
	g(target) = fba.x(target) / 1000;  % specify 0.1% std

	sol = bmfa(model, y, g, sdx, Nsamples, Nchains, sampler, T, svdtol);

end


