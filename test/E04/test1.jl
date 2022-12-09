using LinearAlgebra

function adam(alpha, beta1, beta2, f, grad_f, theta0, eps, convergence_eps)
	m0 = 0
	v0 = 0
	t  = 0
	theta_t = theta0
	theta_prev = theta_t
	mt = m0
	vt = v0

	while true
		t += 1
		gt = grad_f(theta_t)
		mt = beta1 * mt + (1 - beta1) * gt
		vt = beta2 * vt + (1 - beta2) * gt^2
		mt_hat = mt / (1 - beta1^t)
		vt_hat = vt / (1 - beta2^t)
		theta_prev = theta_t
		theta_t = theta_t - (alpha * mt_hat / (sqrt(vt_hat) + eps))

		norm(theta_t - theta_prev) >= convergence_eps || break
	end

	println("Took $t iterations")

	return theta_t
end

f(x) = x^2 - 2x - 1
grad_f(x) = 2x - 2
theta0 = 0

println(adam(0.3, 0.9, 0.999, f, grad_f, theta0, 10e-8, 10e-8))
