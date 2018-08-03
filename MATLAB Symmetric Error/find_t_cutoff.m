function t_star = find_t_cutoff(normhatphiW, tt)
    ind = find(tt >= 0);
    d_phi_W = normhatphiW(ind(2:end)) - normhatphiW(ind(1:end-1));

    if length(find(d_phi_W >= 0)) == 0
        t_star = tt(end);
    else
        first_min_ind = ind(min(find(d_phi_W >= 0)));
        phi_W_threshold = max(normhatphiW(ind(ind >= first_min_ind)));
        tmp = tt(normhatphiW <= phi_W_threshold);
        t_star = min(tmp(tmp>0));
    end
end