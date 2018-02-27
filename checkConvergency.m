function converged = checkConvergency(residuals)
	converged = norm(residuals)<1e-7;
end
	