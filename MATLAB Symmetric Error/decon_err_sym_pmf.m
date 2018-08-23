%% Deconvolution of data X

function [Q, tt, normhatphiW] = decon_err_sym_pmf(W, m, n_tp_iter, n_var_iter)
    n = length(W);
    
    % Precalculate phi_W -------------------------------------------------------
    length_tt = 100;
    mu_K2 = 6;
    RK = 1024 / 3003 / pi;
    hnaive = ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(var(W)) * n^(-1/5);
    hmin = hnaive/3;
    tt = linspace(-1/hmin, 1/hmin, length_tt);

    [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    normhatphiW = sqrt(rehatphiW.^2 + imhatphiW.^2)';
    t_star = find_t_cutoff(normhatphiW, tt);
    tt = linspace(-t_star, t_star, length_tt);

    [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    hat_phi_W = complex(rehatphiW, imhatphiW).';
    % sqrt_psi_hat_W = sqrt(rehatphiW.^2 + imhatphiW.^2)';
    [~, ~, sqrt_psi_hat_W] = compute_psi_W(tt, W);
    sqrt_psi_hat_W = sqrt_psi_hat_W';

    weight_type = 'Epanechnikov';
    weight = KernelWeight(weight_type,tt);

    %--------------------------------------------------------------------------%
    % Solve optimization problem to find PMF
    %--------------------------------------------------------------------------%

    options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
    %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-5);
    %options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);
    
    [A, B] = create_bound_matrices(W, m);

    % ------------------
    % Min T(p)
    % ------------------
    disp(join(["Minimizing T(p)"]))
    drawnow

    fmax = Inf;
    counter = 0;
    while counter < n_tp_iter
        pjtest = unifrnd(0, 1, [1, m]);
        pjtest = pjtest / sum(pjtest);
        xjtest = sort(unifrnd(min(W), max(W), [1, m]));

        [pjnew, xjnew, fmaxnew, exitflag] = min_phase(m, W, pjtest, xjtest, tt, hat_phi_W, ...
                                                      sqrt_psi_hat_W, weight);

        disp(num2str(exitflag))

        if fmaxnew < fmax && exitflag >= 0
            fmax = fmaxnew;
            pj = pjnew;
            xj = xjnew;
            % tp = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
            % [penalty1_max, penalty2_max, ~] = penalties(pj, xj, tt, hat_phi_W);
            disp(join(["T(p) objective =", num2str(fmax)]))
        end

        if exitflag < 0
            disp(join(["Fail with exitflag", num2str(exitflag)]))
        end

        counter = counter + 1;

        drawnow
    end

    % ------------------
    % Min Var
    % ------------------
    
    %Calculate penalties once 
    tp_max = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    [penalty1_max, penalty2_max, ~] = penalties(pj, xj, tt, hat_phi_W);
    disp(join(["tp_max =", num2str(tp_max)]))
    disp(join(["penalties =", num2str(penalty1_max), ",", num2str(penalty2_max)]))

    %Find initial value for varmin based on best solution for fmax
    varmin = var_objective([pj(1:end-1), xj]');
    varmin_init = varmin;

    func = @(x) var_objective(x);
    nonlcon = @(x)phaseconstraint(x, tp_max,penalty1_max, penalty2_max,tt,hat_phi_W,sqrt_psi_hat_W,weight);

    counter = 0;
    
    disp("Minimizing Variance")
    while counter < n_var_iter
        pj_0 = unifrnd(0,1,[1,m]);
        pj_0 = pj_0 / sum(pj_0);
        xj_0 = sort(unifrnd(min(W), max(W), [1,m]));

        x0 = [pj_0(1:end-1), xj_0]';
        [x, fval, exitflag] = fmincon(func, x0, A, B, [], [], [], [], nonlcon, options);
        [pj_new, xj_new] = x_to_pmf(x);

        disp(num2str(exitflag))

        if fval < varmin && exitflag >= 0
            varmin = fval;
            pj = pj_new;
            xj = xj_new;
            disp(num2str(fval))
        end

        if (exitflag >= 0)
            counter = counter + 1;
        end

        drawnow
    end

    tp_final = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    [penalty1_final, penalty2_final, ~] = penalties(pj, xj, tt, hat_phi_W);
    var_final = var_objective([pj(1:end-1), xj]');

    disp(join(["Initial variance was", num2str(varmin_init)]))
    disp(join(["Final variance is", num2str(var_final)]))
    disp(join(["T(p) =", num2str(tp_final)]))
    disp(join(["penalties =", num2str(penalty1_final), ",", num2str(penalty2_final)]))

    % Finalize -----------------------------------------------------------------
    % [pj,xj] = simplify_masses(pj,xj);
    Q.Support = xj;
    Q.ProbWeights = pj;
end

%-------------------------------------------------------------------------------
% Local Functions
%-------------------------------------------------------------------------------

function [pj, xj] = x_to_pmf(x)
    m = (length(x) + 1) / 2;
    x = x';
    pj = [x(1:m-1), 1 - sum( x(1:m-1))];
    xj = x(m:(2 * m - 1));
end

%------------------------%

function [pj, xj, fval, exitflag] = min_phase(m, W, pj, xj, t, hat_phi_W, sqrt_psi_hat_W, weight)
    x = [pj(1:end-1), xj]';
    options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
    %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6);
    % options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);

    [A, B] = create_bound_matrices(W, m);

    func = @(x)tp_objective(x,m,t,hat_phi_W,sqrt_psi_hat_W,weight);
    [x,fval,exitflag] = fmincon(func, x, A, B,[],[],[],[],[],options);

    [pj, xj] = x_to_pmf(x);
end

function [fval, penalty1, penalty2, tp] = tp_objective(x, m, tt, hat_phi_W, sqrt_psi_hat_W, weight)
    %Extract weights and roots
    [pj, xj] = x_to_pmf(x);

    %Integrate
    tp = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);

    %Add penalty terms
    [penalty1, penalty2, ~] = penalties(pj,xj,tt,hat_phi_W);
    penalty_scale = 500;
    fval = tp + penalty_scale * (penalty1 + penalty2);
end

%------------------------%

% function [pj, xj, fval, exitflag] = min_var(m, W, pj, xj, tp_max, penalty1_max, penalty2_max, t, hat_phi_W, sqrt_psi_hat_W, weight)
%     x = [pj(1:end-1),xj]';
%     options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
%     %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-5);
%     %options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);

%     [A, B] = create_bound_matrices(W, m);

%     func = @(x) var_objective(x);
%     nonlcon = @(x)phaseconstraint(x,m,tp_max,penalty1_max, penalty2_max,t,hat_phi_W,sqrt_psi_hat_W,weight);
%     [x,fval,exitflag] = fmincon(func, x, A, B, [], [], [], [], nonlcon, options);

%     [pj, xj] = x_to_pmf(x);
% end

function var1 = var_objective(x)
    [pj, xj] = x_to_pmf(x);
    mean1 = sum(pj * xj');
    var1 = sum(pj * (xj - mean1)'.^2);
end

function [c,ceq] = phaseconstraint(x, tp_max, penalty1_max, penalty2_max, t, hat_phi_W, sqrt_psi_hat_W, weight)
    [pj, xj] = x_to_pmf(x);
    tp = calculate_tp(t,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    [penalty1, penalty2, ~] = penalties(pj,xj,t,hat_phi_W);

    c = [tp - tp_max, penalty1 - penalty1_max, penalty2 - penalty2_max];
    % c = [tp - tp_max, penalty1];
    ceq=[];
end

%------------------------%

function tp = calculate_tp(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight)

    %Calculate characteristic function of our discrete distribution
    [re_phi_p, im_phi_p, norm_phi_p] = computephiX(t, xj, pj);
    phi_p = complex(re_phi_p, im_phi_p);

    %Calculate integrand
    integrand = abs(hat_phi_W - sqrt_psi_hat_W .* phi_p' ./ norm_phi_p').^2.*weight;
    dt = t(2) - t(1);
    tp = dt * sum(integrand);
end

function [penalty1, penalty2, penalty3] = penalties(pj, xj, t, hat_phi_W)

re_hat_phi_W = real(hat_phi_W);
im_hat_phi_W = imag(hat_phi_W);
norm_hat_phi_W = sqrt(re_hat_phi_W.^2 + im_hat_phi_W.^2);
[re_phi_p, im_phi_p, norm_phi_p] = computephiX(t, xj, pj);

%Need phi_U to be real
a = re_phi_p' .* im_hat_phi_W;
b = im_phi_p' .* re_hat_phi_W;
c = a - b;
penalty1  = sum(abs(c));
% penalty1 = 0;

%impose a penalty if |phi_U| is greater than 1:
hat_phi_U = norm_hat_phi_W ./ norm_phi_p';
penalty2 = sum(hat_phi_U(hat_phi_U > 1));
% penalty2 = 0;

penalty3 = 0;   %Removed this penalty
end

function [pj, xj] = simplify_masses(pj, xj)
    zero_tolerance = 0.001;
    search_size = 0.001;

    %Get rid of zeros
    index = 1;
    looping = true;
    while looping
        if pj(index) < zero_tolerance
            xj(index:end-1) = xj(index+1:end);
            xj(end) = [];
            pj(index:end-1) = pj(index+1:end);
            pj(end) = [];
        else
            index = index+1;
        end
        if index >= length(xj)
            looping = false;
        end
    end

    if pj(end) < zero_tolerance
        pj(end) = [];
        xj(end) = [];
    end

    %Get rid of duplicates
    index = 1;
    looping = true;
    if length(xj) < 2
        looping = false;
    end
    while looping
        if abs(xj(index) - xj(index+1)) < search_size
            xj(index) = (xj(index)*pj(index)+xj(index+1)*pj(index+1))/(pj(index)+pj(index+1));  %Weighted average location
            xj(index+1:end-1) = xj(index+2:end);
            pj(index) = pj(index)+pj(index+1);
            pj(index+1:end-1) = pj(index+2:end);
            pj(end) = [];
            xj(end) = [];
        else
            index = index+1;
        end
        if index == length(xj)
            looping = false;
        end
    end

    pj = pj/sum(pj);    %Normalise
end

function weight = KernelWeight(weight_type, x)
    length_x = length(x);
    switch weight_type
        case 'Epanechnikov'
            sig_x=-x(1)/2;
            weight=0.75/(2*sig_x)*(1-(x/(2*sig_x)).^2);
        case 'Uniform'
            max_x=-x(1);
            weight = zeros(1,length_x) + 1/(2*max_x);
        case 'Triangular'
            max_x = -x(1);
            weight = -1/(max_x^2)*(abs(x)-max_x);
        case 'Triweight'
            max_x = -x(1);
            weight = 35/(32*max_x)*(1 - (x/max_x).^2).^3;
    end
end

function [A, B] = create_bound_matrices(W, m)
    % pj non-negative
    A = zeros(2*m + 1, 2*m - 1);
    A(1:(m-1), 1:(m-1)) = eye( m - 1 );
    B = zeros(2 * m - 1, 1);
    % pj sum to less than 1
    A(m, 1:m-1) = -ones(1,m-1);
    B(m) = -1;

    % thetaj are increasing
    for i = 1:(m-1)
        A(m+i, (m+i-1):(m+i)) = [-1, 1];
    end
    % min(W) < thetaj < max(W)
    A(2*m, m) = 1;
    A(2*m+1, 2*m-1) = -1;
    B(2*m) = min(W);
    B(2*m+1) = -max(W);

    A = -A;
    B = -B;
end