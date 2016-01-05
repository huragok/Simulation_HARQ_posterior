function E = get_factor_PEP_update_param(dist_sqr, beta, K, sigma_sqr, x)

tmp0 = x * sigma_sqr * (K + 1);
E = tmp0 ./ (tmp0 + beta * dist_sqr) .* exp(-(K * beta * dist_sqr) ./ (tmp0 + beta * dist_sqr)); % Note the difference between the exponential integral function definition in the reference and Matlab

