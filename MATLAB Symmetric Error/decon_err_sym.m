function [fXdeconvoluted, xx, Q] = decon_err_sym(W, xx, m, show_diagnostics)

    if (~exist('xx','var') | isempty(xx))
        xx = linspace(min(W), max(W), 100);
    end

    if (~exist('m', 'var') | isempty(m))
      m = 20;
    end

    if (~exist('show_diagnostics', 'var') | isempty(show_diagnostics))
      show_diagnostics = 0;
    end

    % Deconvolve to pmf --------------------------------------------------------
    n_tp_iter = 5;
    n_var_iter = 5;
    [Q, tt, normhatphiW] = decon_err_sym_pmf(W, ...
                                             m, ...
                                             n_tp_iter, ...
                                             n_var_iter,...
                                             show_diagnostics);

    % Convert pmf to pdf -------------------------------------------------------    
    fXdeconvoluted = decon_err_sym_pmf2pdf(xx, ...
                                           tt, ...
                                           Q.Support, ...
                                           Q.ProbWeights, ...
                                           W);
end