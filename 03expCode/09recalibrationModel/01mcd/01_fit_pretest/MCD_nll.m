function nll = MCD_nll(s_unique_sec, duration, fs, onset, stim_dura, ...% fixed parameters
    nT_vfirst, nT_simul, nT_afirst,...% behavioral responses to be fitted
    p)% free parameters

[P_AFIRST, P_VFIRST, P_SIMUL] = MCD_prob(-s_unique_sec, duration, fs, onset, stim_dura, ...
    p(1), p(2), p(3), p(4), p(5), p(6));

nll = -nT_afirst * log(P_AFIRST)' - nT_vfirst * log(P_VFIRST)' - nT_simul * log(P_SIMUL)';

end
