function X = conditional_mvn(X, mu, sigma, inds)
% extrapolates from the 'inds' rows of X~N(mu,sigma) the remaining rows	
% i.e. distribution of unbounded (u) given bounded (b)

	I = inds;
	Cub = sigma(~I,I);
	Cbu = sigma(I,~I);
	Cbb = sigma(I,I);
	Cuu = sigma(~I,~I);
	B = Cub/Cbb;
	fcov = Cuu - B*Cbu;
	fcov = 0.5*(fcov + fcov'); % make symmetric
	mu_b = mu(I);
	mu_u = mu(~I);
	
	% loop version
%	for i=1:size(X,2)
%		for c=1:size(X,3)
%			v = X(I,i,c);
%			X(~I,i,c) = mvnrnd(mu_u + B*(v-mu_b), fcov)'; 
%		end
%	end
	
	% compute all conditionals at the same time (shared covariance)
	% we want to compute: mvnrnd(mu_u + B*(v-mu_b), fcov)
	Nd = length(mu_u);
	Ns = size(X,2);
	Nc = size(X,3);
	Z = mvnrnd(mu_u, fcov, Nc*Ns)';
	Z = reshape(Z, Nd, Ns, Nc);
	
	% add the missing mean term
	for i=1:Ns
		for c=1:Nc
			v = X(I,i,c);
			X(~I,i,c) = Z(:,i,c) + B*(v-mu_b);
		end
	end
	
end

