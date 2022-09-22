function p_tilde = add_lapse(p1, p2, p3, epsilon, kappa)
    p_tilde = (1-epsilon) * p1 + epsilon * kappa * p2 + epsilon * kappa * p3;
end