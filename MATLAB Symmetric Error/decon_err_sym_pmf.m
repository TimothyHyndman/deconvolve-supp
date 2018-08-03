%% Deconvolution of data X

function [Q, tt, normhatphiW] = decon_err_sym_pmf(W, m, n_tp_iter, n_var_iter)
    n = length(W);
    
    % Precalculate phi_W -------------------------------------------------------
    length_tt = 100;
    mu_K2 = 6;
    RK = 1024 / 3003 / pi;
    hnaive = ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * var(W) * n^(-1/5);
    hmin = hnaive/3;
    tt = linspace(-1/hmin, 1/hmin, length_tt);

    [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    normhatphiW = sqrt(rehatphiW.^2 + imhatphiW.^2)';
    t_star = find_t_cutoff(normhatphiW, tt);




    % %Refined interval of t-values 
    % tt = linspace(-t_star, t_star, length_tt+1); 
    % [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    % hat_phi_W = complex(rehatphiW, imhatphiW)';

    % [rehatpsi, imhatpsi, sqabshatpsi] = compute_psi_W(tt, W);
    % % sqrt_psi_hat_W = abs(hat_phi_W);
    % sqrt_psi_hat_W = sqabshatpsi';


    tt1 = -t_star; 
    tt2 = t_star; 
     tt=tt1:(tt2-tt1)/length_tt:tt2; %Refined interval of t-values  
    n = length(W); 
    hat_phi_W = 0; 
    for i=1:n 
        hat_phi_W = hat_phi_W + exp(1i*tt*W(i)); 
    end 
    hat_phi_W = (1/n)*hat_phi_W; 
    sqrt_psi_hat_W = abs(hat_phi_W); 






    % Minimize difference in phase functions -----------------------------------
    weight_type = 'Epanechnikov';
    weight = KernelWeight(weight_type,tt);

    fmax = Inf;
    n_iterations = n_tp_iter;
    n_var_iterations = n_var_iter;
    
    for i = 1:n_iterations
        pjtest = unifrnd(0, 1, [1, m]);
        pjtest = pjtest / sum(pjtest);
        xjtest = sort(unifrnd(min(W), max(W), [1, m]));

        [pjnew,xjnew,fmaxnew,exitflag] = min_phase(m,W,pjtest,xjtest,tt,hat_phi_W,sqrt_psi_hat_W,weight);

        if fmaxnew < fmax && exitflag>0
            fmax = fmaxnew;
            pj = pjnew;
            xj = xjnew;
            display("pass")
        end

        if exitflag <= 0
            display("fail")
        end
    end

    % Minimize variance --------------------------------------------------------
    %Calculate penalties once 
    tp_max = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight)
    [penalty1_max, penalty2_max, ~] = penalties(pj, xj, tt, hat_phi_W);

    %Find initial value for varmin based on best solution for fmax
    varmin = var_objective([pj(1:end-1),xj]')

    for i = 1:n_var_iterations
        pjtest = unifrnd(0,1,[1,m]);
        pjtest = pjtest / sum(pjtest);
        xjtest=sort(unifrnd(min(W), max(W), [1,m]));
        [pjnew,xjnew,fval,exitflag] = min_var(m,W,pjtest,xjtest,tp_max,penalty1_max,penalty2_max,tt,hat_phi_W,sqrt_psi_hat_W,weight);
        if fval < varmin && exitflag > 0
            varmin = fval;
            pj = pjnew;
            xj = xjnew;
        end   
    end

    tp_max = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight)
    varmin
    [penalty1_max, penalty2_max, ~] = penalties(pj, xj, tt, hat_phi_W);
    

    % fval2 = tp_objective([pj,xj],m,tt,hat_phi_W,sqrt_psi_hat_W,weight);
    % if fval2/tp_max < 0.999
    %     display('T(Y) got smaller!')
    % end

    % Finalize -----------------------------------------------------------------
    [pj,xj] = simplify_masses(pj,xj);
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

%------------------------%

function [pj, xj, fval, exitflag] = min_var(m, W, pj, xj, tp_max, penalty1_max, penalty2_max, t, hat_phi_W, sqrt_psi_hat_W, weight)
    x = [pj(1:end-1),xj]';
    options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
    %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-5);
    %options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);

    % pj non-negative
    A_tp = zeros(2*m + 1, 2*m - 1);
    A_tp(1:(m-1), 1:(m-1)) = eye( m - 1 );
    B_tp = zeros(2 * m - 1, 1);
    % pj sum to less than 1
    A_tp(m, 1:m-1) = -ones(1,m-1);
    B_tp(m) = -1;
    % thetaj are increasing
    for i = 1:(m-1)
        A_tp(m+i, (m+i-1):(m+i)) = [-1, 1];
    end
    % min(W) < thetaj < max(W)
    A_tp(2*m, m) = 1;
    A_tp(2*m+1, 2*m-1) = -1;
    B_tp(2*m) = min(W);
    B_tp(2*m+1) = -max(W);

    A_tp = -A_tp;
    B_tp = -B_tp;

    func = @(x) var_objective(x);
    nonlcon = @(x)phaseconstraint(x,m,tp_max,penalty1_max, penalty2_max,t,hat_phi_W,sqrt_psi_hat_W,weight);
    [x,fval,exitflag] = fmincon(func,x,A_tp,B_tp,[],[],[],[],nonlcon,options);

    [pj, xj] = x_to_pmf(x);
end

function fval = var_objective(x)
    [pj, xj] = x_to_pmf(x);
    fval = sum(pj.*(xj.^2)) - (sum(pj.*xj))^2;
end

function [c,ceq] = phaseconstraint(x, m, tp_max, penalty1_max, penalty2_max, t, hat_phi_W, sqrt_psi_hat_W, weight)
    [pj, xj] = x_to_pmf(x);

    %My own integral thingo
    dt = t(2) - t(1);
    tp = calculate_tp(t,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);

    [penalty1, penalty2, ~] = penalties(pj,xj,t,hat_phi_W);

    %Finish off
    tol = 0;
    c = [tp - tp_max; penalty1 - penalty1_max; penalty2 - penalty2_max] - tol;
    % c = tp - tp_max - tol;
    ceq=[];
end

%------------------------%

function [pj, xj, fval, exitflag] = min_phase(m, W, pj, xj, t, hat_phi_W, sqrt_psi_hat_W, weight)
    x = [pj(1:end-1), xj]';
    options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
    %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6);
    %options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);

    % pj non-negative
    A_tp = zeros(2*m + 1, 2*m - 1);
    A_tp(1:(m-1), 1:(m-1)) = eye( m - 1 );
    B_tp = zeros(2 * m - 1, 1);
    % pj sum to less than 1
    A_tp(m, 1:m-1) = -ones(1,m-1);
    B_tp(m) = -1;
    % thetaj are increasing
    for i = 1:(m-1)
        A_tp(m+i, (m+i-1):(m+i)) = [-1, 1];
    end
    % min(W) < thetaj < max(W)
    A_tp(2*m, m) = 1;
    A_tp(2*m+1, 2*m-1) = -1;
    B_tp(2*m) = min(W);
    B_tp(2*m+1) = -max(W);

    A_tp = -A_tp;
    B_tp = -B_tp;

    func = @(x)tp_objective(x,m,t,hat_phi_W,sqrt_psi_hat_W,weight);
    [x,fval,exitflag] = fmincon(func, x, A_tp, B_tp,[],[],[],[],[],options);

    [pj, xj] = x_to_pmf(x);
end

function [fval, penalty1, penalty2, tp] = tp_objective(x, m, tt, hat_phi_W, sqrt_psi_hat_W, weight)
    %Extract weights and roots
    [pj, xj] = x_to_pmf(x);

    %Integrate
    dt = tt(2) - tt(1);
    tp = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);

    %Add penalty terms
    [penalty1, penalty2, ~] = penalties(pj,xj,tt,hat_phi_W);
    fval = tp + penalty1 + penalty2;
end

function tp = calculate_tp(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight)

    %Calculate characteristic function of our discrete distribution
    m = length(pj);
    phi_p = 0;
    for i=1:m
        phi_p = phi_p + (pj(i)*exp(1i*t*xj(i)));
    end

    %Calculate integrand
    integrand = abs(hat_phi_W.*conj(phi_p) - abs(phi_p).*sqrt_psi_hat_W).^2.*weight;
    dt = t(2) - t(1);
    tp = dt * sum(integrand);
end

function [penalty1, penalty2, penalty3] = penalties(pj, xj, t, hat_phi_W)

%Calculate characteristic function of our discrete distribution
m = length(pj);
phi_p = 0;
for i=1:m
    phi_p = phi_p + (pj(i)*exp(1i*t*xj(i)));
end
%phi_p = conj(pj*(exp(1i*t'*xj)')); %Slightly slower

re_phi_p = real(phi_p);
im_phi_p = imag(phi_p);
re_hat_phi_W = real(hat_phi_W);
im_hat_phi_W = imag(hat_phi_W);
norm_hat_phi_W = sqrt(re_hat_phi_W.^2+im_hat_phi_W.^2);
norm_phi_p = sqrt(re_phi_p.^2+im_phi_p.^2);

%Need phi_U to be real
penalty1=sum(abs(re_phi_p .* im_hat_phi_W - im_phi_p .* re_hat_phi_W));

%Want phi_W bar(phi_p) = |phi_W||phi_p|
% error_vec = hat_phi_W.*conj(phi_p) - norm(hat_phi_W).*norm(phi_p);
% %error_vec = hat_phi_W - norm(hat_phi_W).*phi_p/(norm(phi_p));
% penalty1 = sum(error_vec);

%impose a penalty if |phi_U| is greater than 1:
tolerance = 0;
hat_phi_U=norm_hat_phi_W ./ (norm_phi_p);
penalty2=sum(hat_phi_U(hat_phi_U > 1 + tolerance));

%impose a penalty if phi_U is < 0
re_phi_U = real(hat_phi_U);
penalty3 = -sum(re_phi_U(re_phi_U < 0 - tolerance));

%Scale
scale = 1;
penalty1 = scale*penalty1;
penalty2 = scale*penalty2;
end

%------------------------%

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
