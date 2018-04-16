function vmap = tmvn_map(mu, sigma, lb, ub, mode)
% finds maximal value of truncated MVN using quadprog
	
	disp('Finding initial flux...');
	d = length(mu);

	if ~exist('mode','var') % 'inv' or 'chol', how to compute the mode
		mode = 'chol';
	end

	if ~exist('ub','var')
		ub = inf(d,1);
	end
	if ~exist('lb','var')
		lb = -inf(d,1);
	end

	% centerize
	lb = lb - mu;
	ub = ub - mu;
	init = max(lb,0) + 0.01;
	tol = 1e-8;
	options = optimoptions('quadprog', 'Display', 'none', 'TolX', tol, 'TolFun', tol, 'TolCon', tol);
	
	switch mode
		case 'inv'
			iS = inv(sigma);
			iS = 0.5*(iS + iS');
			
			vmap = quadprog(iS, [], [], [], [], [], lb, ub, init, options);
			vmap = vmap + mu;
		case 'chol'
			% whitened domain
			L = chol(sigma,'lower');
			g = [-lb; ub];
			trueb = ~isinf(g);   % remove NaN's 
			g = g(trueb);

			F = [eye(d); -eye(d)];
			F = F(trueb,:);
			F = F*L;

			vmap = quadprog(eye(size(L,2)), [], -F, g, [], [], [], [], L\init, options);
			vmap = L*vmap + mu;
	end
end

