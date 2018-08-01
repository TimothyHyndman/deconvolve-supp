%% Deconvolution of data X

function [Q,tt,tt1,tt2,normhatphiW] = decon_err_sym_pmf(W)
    %Basic parameters
    n = length(W);
    m = 10;

    %-------------------------------------------------------------
    %-------------------------------------------------------------
    % Precalculate phi_W
    %-------------------------------------------------------------
    %-------------------------------------------------------------

    length_tt=100;
    a=-8;
    b=8;
    tt=a:((b-a)/length_tt):b;
    [tt,tt1,tt2,hat_phi_W,sqrt_psi_hat_W,normhatphiW]=computephiW(tt,length_tt,W,n);

    %Choose kernel Weight
    weight_type = 'Epanechnikov';
    weight = KernelWeight(weight_type,tt);

    %-------------------------------------------------------------
    %-------------------------------------------------------------
    % Solve for xj,pj
    %-------------------------------------------------------------
    %-------------------------------------------------------------

    fmax = 10^6;   %Just a big number
    n_iterations = 5; %Setting this to 1 works perfectly almost all the time.
    n_var_iterations = 5;  %Usually gets there in first iteration but sometimes needs a few goes

    %Minimize difference in phase functions
    for i = 1:n_iterations
        pjtest = unifrnd(0,1,[1,m]);
        pjtest = pjtest/sum(pjtest);
        xjtest=sort(unifrnd(min(W),max(W),[1,m]));
        [pjnew,xjnew,fmaxnew,exitflag] = min_phase(m,W,pjtest,xjtest,tt,hat_phi_W,sqrt_psi_hat_W,weight);
        %If our new answer is better than our old one then accept it
        if fmaxnew < fmax && exitflag>0
            fmax = fmaxnew;
            pj = pjnew;
            xj = xjnew;
        end
    end

    %Calculate penalties once 
    dt = tt(2) - tt(1);
    integrand = insideintegral(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    fmax0 = dt*sum(integrand);

    %Find initial value for varmin based on best solution for fmax
    varmin = var_objective([pj,xj]);

    for i = 1:n_var_iterations
        pjtest = unifrnd(0,1,[1,m]);
        pjtest = pjtest/sum(pjtest);
        xjtest=sort(unifrnd(min(W),max(W),[1,m]));
        [pjnew,xjnew,fval,exitflag] = min_var(m,W,pjtest,xjtest,fmax0,tt,hat_phi_W,sqrt_psi_hat_W,weight);
        if fval < varmin && exitflag > 0
            varmin = fval;
            pj = pjnew;
            xj = xjnew;
        end   
    end

    fval2 = phase_objective([pj,xj],m,tt,hat_phi_W,sqrt_psi_hat_W,weight);
    if fval2/fmax0 < 0.999
        display('T(Y) got smaller!')
    end

    [pj,xj] = simplify_masses(pj,xj);
    Q.Support = xj;
    Q.ProbWeights = pj;
end

%-------------------------------------------------------------------------------
% Local Functions
%-------------------------------------------------------------------------------

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

function [tt, tt1, tt2, hat_phi_W, sqrt_psi_hat_W, normhatphiW] = computephiW(tt, length_tt, W, n)
    OO=outerop(tt,W,'*');

    %-----------------------------------------------
    %Estimate empirical characersitic fucntion of W
    %-----------------------------------------------

    rehatphiW=sum(cos(OO),2)/n;
    imhatphiW=sum(sin(OO),2)/n;
    normhatphiW=sqrt(rehatphiW.^2+imhatphiW.^2);

    %------------------------------------------------------------
    %t^*: Keep only t-values such that |\hat_phi_W|< n^{-0.25}
    %------------------------------------------------------------

    tmp=tt(normhatphiW<n^(-0.25));
    %Refine the interval of t-values accordingly
    tt1=max(tmp(tmp<0));
    tt2=min(tmp(tmp>0));
    if isempty(tmp(tmp<0))
        tt1=min(tt);
    end

    if isempty(tmp(tmp>0))
        tt2=max(tt);
    end
    
    tt=tt1:(tt2-tt1)/length_tt:tt2; %Refined interval of t-values 
    hat_phi_W = 0;
    for i=1:n
        hat_phi_W = hat_phi_W + exp(1i*tt*W(i));
    end
    hat_phi_W = (1/n)*hat_phi_W;
    sqrt_psi_hat_W = abs(hat_phi_W);
end

function y = outerop(a, b, operator)
    if nargin<3
        operator='+';                       % for only two arguments assume outerproduct 
    end  

    if isequal(operator,'*')                % common operator 
        y=a(:)*b(:)';
    else    
      outera=a(:)*ones(1,length(b));        % these are the matrices that 
      outerb=ones(length(a),1)*b(:).';      % meshgrid(A,B) would generate 
      functionHandle=str2func(operator);  % new R14 functionality
      y=functionHandle(outera,outerb);    % allows faster/neater method
    end
end

%------------------------%

function [pj, xj, fval, exitflag] = min_var(m, W, pj, xj, fmax0, t, hat_phi_W, sqrt_psi_hat_W, weight)
    x = [pj,xj];
    options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
    %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-5);
    %options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);

    %Setup constraints in form A*x<b. We need that entries of pj are
    %non-negative
    A = zeros(2*m);
    A(1:m,1:m) = -eye(m);
    b = zeros(2*m,1);

    %Optional: The xj are increasing
    for i = 1:m-1
        A(m+i,m+i) = 1;
    end
    for i = 1:m-1
        A(m+i,m+i+1) = -1;
    end
    b(2*m) = 1;

    %Setup constraints in the form Aeq*x = beq. We need that the pj add up to one.
    Aeq = zeros(2*m);
    beq = zeros(2*m,1);
    Aeq(2*m,1:m) = ones(1,m);
    beq(2*m) = 1;

    %Setup bounds lb<x<ub
    %Masses
    lb = zeros(2*m,1);
    ub = ones(2*m,1);
    %Roots
    lb(m+1:2*m) = min(W);
    ub(m+1:2*m) = max(W);
    func = @(x) var_objective(x);
    nonlcon = @(x)phaseconstraint(x,m,fmax0,t,hat_phi_W,sqrt_psi_hat_W,weight);
    [x,fval,exitflag] = fmincon(func,x,A,b,Aeq,beq,lb,ub,nonlcon,options);

    pj = x(1:m);
    xj = x(m+1:2*m);
end

function fval = var_objective(x)
    m = length(x)/2;
    pj = x(1:m);
    xj = x(m+1:2*m);
    fval = sum(pj.*(xj.^2)) - (sum(pj.*xj))^2;
end

function [c,ceq] = phaseconstraint(x, m, fmax0, t, hat_phi_W, sqrt_psi_hat_W, weight)
    pj = x(1:m);
    xj = x(m+1:2*m);

    %My own integral thingo
    dt = t(2) - t(1);
    integrand = insideintegral(t,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    Tp = dt*sum(integrand);

    %Finish off
    tol = 0;
    %c=[Tp - fmax0 - tol; penalty1-fmax1 - tol ; penalty2 - fmax2-tol];
    c = Tp - fmax0 - tol;
    ceq=[];
end

%------------------------%

function [pj, xj, fval, exitflag] = min_phase(m, W, pj, xj, t, hat_phi_W, sqrt_psi_hat_W, weight)
    x = [pj,xj];
    options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6); 
    %options = optimoptions('fmincon','Display','off','Algorithm','interior-point','TolFun',1e-6);
    %options = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-6);

    %Setup constraints in form A*x<b. We need that entries of pj are
    %non-negative
    A = zeros(2*m);
    A(1:m,1:m) = -eye(m);
    b = zeros(2*m,1);

    %Optional: The xj are increasing
    for i = 1:m-1
        A(m+i,m+i) = 1;
    end
    for i = 1:m-1
        A(m+i,m+i+1) = -1;
    end
    b(2*m) = 1;

    %Setup constraints in the form Aeq*x = beq. We need that the pj add up to one.
    Aeq = zeros(2*m);
    beq = zeros(2*m,1);
    Aeq(2*m,1:m) = ones(1,m);
    beq(2*m) = 1;

    %Setup bounds lb<x<ub
    lb = zeros(2*m,1);
    lb(m+1:2*m) = min(W);
    ub = ones(2*m,1);
    ub(m+1:2*m) = max(W);

    func = @(x)phase_objective(x,m,t,hat_phi_W,sqrt_psi_hat_W,weight);
    [x,fval,exitflag] = fmincon(func,x,A,b,Aeq,beq,lb,ub,[],options);

    pj = x(1:m);
    xj = x(m+1:2*m);
end

function [fval, penalty1, penalty2, Tp] = phase_objective(x, m, tt, hat_phi_W, sqrt_psi_hat_W, weight)
    %Extract weights and roots
    pj = x(1:m);
    xj = x(m+1:2*m);

    %Integrate
    dt = tt(2) - tt(1);
    integrand = insideintegral(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    Tp = dt*sum(integrand);

    %Add penalty terms
    % [penalty1,penalty2] = penalties(pj,xj,tt,hat_phi_W);  %My code is working
    % without these!! :)
    penalty1 = 0;
    penalty2 = 0;
    fval = Tp+penalty1+penalty2;
end

function integrand = insideintegral(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight)

    %Calculate characteristic function of our discrete distribution
    m = length(pj);
    phi_p = 0;
    for i=1:m
        phi_p = phi_p + (pj(i)*exp(1i*t*xj(i)));
    end
    %phi_p = conj(pj*(exp(1i*t'*xj)')); %Slightly Slower

    %Calculate integrand
    %fred = hat_phi_W - sqrt_psi_hat_W.*phi_p./abs(phi_p);  %Slightly slower
    fred = hat_phi_W.*conj(phi_p) - abs(phi_p).*sqrt_psi_hat_W;

    integrand = abs(fred).^2.*weight;
end

%------------------------%

function [pj, xj] = simplify_masses(pj, xj)
    zero_tolerance = 0.0001;
    search_size = 0.01;

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
