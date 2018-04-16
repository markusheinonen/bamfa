function Xs = tmvn_gibbs2(mu,sigma,lower,upper,x0,N,C,Nskip, svdtol)
% Gibbs sampler for truncated MVN
% code translated from tmvnorm.rtmvnorm() function by Stefan Wilhelm

	D = length(mu);
	
	if ~exist('Nskip','var') % save only every 10'th sample
		Nskip = 10;
	end
	if ~exist('C','var') % 5 chains
		C = 5;
	end
	if ~exist('N','var') % 100 samples
		N = 100;
	end
	if ~exist('x0','var') % initial value
		x0 = min(max(zeros(D,1), lower),upper);
	end
	if ~exist('svdtol','var') % 100 samples
		svdtol = 0;
	end

	% make zero mean
	lb = lower - mu;
	ub = upper - mu;
	x0 = x0 - mu;
	
	% whiten by cholesky
%	fprintf(' Whitening %d x %d covariance matrix..\n',D,D);
	if svdtol == 0
		L = chol(sigma,'lower');
	else
		% not sure if this is a good idea...
		[L,S] = ldl(sigma,'lower');
		d = diag(S);
		d(d<svdtol) = 0;
		S = diag(d);
		L = L*sqrt(S);
		L = L(:,var(L)>0);
		S = size(L,2);
		
%		[U,S] = eig(sigma);
%		L = U*sqrt(S);
	end
	
	% transform initial value into whitened domain
	x0w = L\x0;
	D = length(x0w);
	
	% initial sample
	Xs = zeros(D,N,C);
	Xs(:,1,:) = repmat(x0w, 1, C);
	X = repmat(x0w, 1, C);

	% relax bounds so that x0 is within them
	lb = lb + min(0, L*x0w - lb);
	ub = ub + max(0, L*x0w - ub);
	
	LB = repmat(lb,1,C);
	UB = repmat(ub,1,C);
	
	% precompute L splits
	Lp = L>0;
	Ln = L<0;
	
	fprintf('Sampling %d chains of %d bounded fluxes... [%d sweeps per sample, %d updates per sample]\n', C, D, Nskip, D*Nskip);

	tic;
	startcpu = cputime;
		
	for j=1:N
		for k=1:Nskip
			LX = L*X;
			
			for i=D:-1:1  % sweep from end
				Li = L(:,i);
				Lip = Lp(:,i);
				Lin = Ln(:,i);
				
				% Li & Ghosh 2015
				lis_lower = (LB - LX + Li*X(i,:)) ./ repmat(Li,1,C);
				lis_upper = (UB - LX + Li*X(i,:)) ./ repmat(Li,1,C);
				lp = max(lis_lower(Lip,:), [], 1);
				ln = max(lis_upper(Lin,:), [], 1);
				up = min(lis_upper(Lip,:), [], 1);
				un = min(lis_lower(Lin,:), [], 1);

				ai = max([lp;ln], [], 1);
				bi = min([up;un], [], 1);

				ai = min(ai,bi); % ensure non-negative range

				old = X(i,:);
%				new = rtnorm(ai,bi);   % much slower than trandn1
%				new = trandn1(ai,bi);
				new = trandn(ai,bi)';
				X(i,:) = new;
				LX = LX + Li*(new - old);  % update the L*X
%				display(sprintf('i=%d:  %.2f --> %.2f  [%.2f %.2f]', i, old, new, ai, bi));
			end
		end
		
		if j <= 5
			fprintf(' %d/%d fluxes sampled in %.1fmin (%.1fmin cpu)\n', j, N, toc/60, (cputime-startcpu)/60);
		elseif j <= 100 && mod(j,10)==0
			fprintf(' %d/%d fluxes sampled in %.1fmin (%.1fmin cpu)\n', j, N, toc/60, (cputime-startcpu)/60);
		elseif mod(j,100)==0
			fprintf(' %d/%d fluxes sampled in %.1fmin (%.1fmin cpu)\n', j, N, toc/60, (cputime-startcpu)/60);
		end

		Xs(:,j,:) = X;
	end
	

	% unwhiten and add mean
%	Xs = L*Xs + repmat(mu,1,N);

	Z = cellfun(@(X) L*X, num2cell(Xs, [1 2]), 'un',0);
	Xs = cat(3,Z{:});
	Xs = Xs + repmat(mu,1,N,C);
end

