function E = get_factor_PEP_update_param(dist_sqr, beta, sigma_sqr, x)

E = 1 ./ (1 + dist_sqr * beta / x / sigma_sqr); % Note the difference between the exponential integral function definition in the reference and Matlab

