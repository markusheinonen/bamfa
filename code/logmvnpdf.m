function [val,term1,term2] = logmvnpdf(x,m,S)
% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% x is NxD, mu is 1xD, Sigma is DxD

	logdet = @(A) 2*sum(log(diag(chol(A))));

	D = size(x,2);
	const = -0.5 * D * log(2*pi);

	term1 = diag(-0.5 * (x-repmat(m,1,D))' / S * (x-repmat(m,1,D)));
	
	try
		term2 = const - 0.5 * logdet(S);
	catch
		warning('chol failed: increasing diagonal to fix');
		fix = abs(min(eig(S)));
%		min(eig( S + fix  ))
		term2 = const - 0.5 * logdet(S + (fix*1.01)*eye(size(S)));
	end
	
	val = term1 + term2;
end

