function fX = decon_err_sym_pmf2pdf(xx, tt, theta, p, W)
    
    % Estimate sd_U ------------------------------------------------------------
    tt_BB_length = 200;
    tt_BB = linspace(tt(1), tt(end), tt_BB_length);
    var_U = estimate_var_u(W, tt_BB, theta, p);

    % Estimate PhiX and PhiU ---------------------------------------------------
    [rephip, imphip, normphip] = computephiX(tt, theta, p);
    [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    normhatphiW = sqrt(rehatphiW.^2 + imhatphiW.^2);
    hatphiU = normhatphiW ./ normphip;

    % Adjust estimator of phi_U as recommended in the paper --------------------
    tlim = [min(tt), max(tt)];
    ppphiU = spline(tt, hatphiU);

    % Find Plug-In Bandwidth ---------------------------------------------------
    h = PI_deconvUestth4(W, tlim, ppphiU, var_U);

    % Compute estimator --------------------------------------------------------
    fX = fXKernDec2(xx, h, W, tlim, ppphiU, var_U);

    %Remove negative parts and rescale to integrate to 1
    fX(fX < 0) = 0 * fX(fX < 0);
    dx = xx(2) - xx(1);
    fX = fX / sum(fX) / dx;
end

% ------------------------------------------------------------------------------
% Local Functions
% ------------------------------------------------------------------------------

function var_U = estimate_var_u(W, tt_BB, theta, p)

    [re_phi_W_BB, im_phi_W_BB] = compute_phi_W(tt_BB, W);
    norm_phi_W_BB = sqrt(re_phi_W_BB.^2 + im_phi_W_BB.^2);

    [re_phi_X_BB, im_phi_X_BB, norm_phi_X_BB] = computephiX(tt_BB, theta, p);
    hatphiUBB = norm_phi_W_BB ./ norm_phi_X_BB;

    t_vec = tt_BB(hatphiUBB' >= 0.95);
    phi_U_t_vec = hatphiUBB(hatphiUBB' >= 0.95);

    pp = polyfit(t_vec, phi_U_t_vec', 2);
    var_U = abs(2*pp(1));
end