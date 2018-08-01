%% Plot deconvolution stuff

function [fXdeconvoluted, xx, Q] = decon_err_sym(W, xx)
    if ~exist('xx','var')
        xx = linspace(min(W), max(W), 100);
    end

    % Deconvolve to pmf --------------------------------------------------------
    [Q, tt, tt1, tt2, normhatphiW] = decon_err_sym_pmf(W);
    

    % Convert pmf to pdf -------------------------------------------------------    
    fXdeconvoluted = makesmooth(xx, tt, tt1, tt2, Q.Support, Q.ProbWeights, W, normhatphiW);
end

% ------------------------------------------------------------------------------
% Local Functions
% ------------------------------------------------------------------------------

function fXdecc = makesmooth(xx, tt, tt1, tt2, xgrid, psol, W, normhatphiW)
    n = length(W);
    dx=xx(2)-xx(1);

    %----------------------------
    % Estimate phi_X and phi_U
    %----------------------------

    [rephip,imphip,normphip]=computephiX(tt,xgrid,psol);
    hatphiU=normhatphiW./normphip;

    %---------------------------------------------------------------------------
    % Estimate var(U): approximate phi_U by poly of degree 2, and estimate varU 
    % by -2 phi_U''
    %---------------------------------------------------------------------------

    %For this we use a finer grid than the grid tt
    ttBB=tt1:(tt2-tt1)/200:tt2;
    OO=outerop(ttBB,W,'*');

    %Estimate empirical characersitic function of W
    rehatphiWBB=sum(cos(OO),2)/n;
    imhatphiWBB=sum(sin(OO),2)/n;

    clear OO;
    normhatphiWBB=sqrt(rehatphiWBB.^2+imhatphiWBB.^2);
    clear rehatphiWBB imhatphiWBB;

    [rephipBB,imphipBB,normphipBB]=computephiX(ttBB,xgrid,psol);
    hatphiUBB=normhatphiWBB./normphipBB;

    tvec=ttBB(hatphiUBB'>=1-0.05);
    phiUtvec=hatphiUBB(hatphiUBB'>=1-0.05);
    pp=polyfit(tvec,phiUtvec',2);
    hatvarU=-2*pp(1);

    %Then compute density estimator as indicated in the paper

    tlim=[min(tt),max(tt)];

    %Adjust estimator of phi_U as recommended in the paper
    ppphiU=spline(tt,hatphiU);

    %Compute bandwidth as recommended in the paper

    hPIc=PI_deconvUestth4(W,tlim,ppphiU,hatvarU);

    %Compute density estmator
    fXdecc=fXKernDec2(xx,hPIc,W,tlim,ppphiU,hatvarU);

    %Remove negative parts and rescale to integrate to 1
    fXdecc(fXdecc<0)=0*fXdecc(fXdecc<0);
    fXdecc=fXdecc/sum(fXdecc)/dx;
end

function [rephip, imphip, normphip] = computephiX(tt, xgrid, psol)
    OO=outerop(tt,xgrid,'*');
    pmat=repmat(psol,length(tt),1);
    cosO=cos(OO).*pmat;
    sinO=sin(OO).*pmat;
    clear OO;

    rephip=sum(cosO,2);
    imphip=sum(sinO,2);
    normphip=sqrt(rephip.^2+imphip.^2);
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

function fXdec = fXKernDec2(xx, hPI, W, tlim, ppphiU, hatvarU)
	phiK = @(t) (1-t.^2).^3;
	muK2 = 6;

	dt = .0002;
	t = (-1:dt:1)';
	t=reshape(t,length(t),1);
	
	phiU=phiUspline(t/hPI,hatvarU,tlim,ppphiU);
	OO=outerop(t/hPI,W,'*');
	%Estimate empirical characersitic fucntion of W
	
	n=length(W);
	rehatphiW=sum(cos(OO),2)/n;
	imhatphiW=sum(sin(OO),2)/n;

	rephip=rehatphiW./phiU;
	imphip=imhatphiW./phiU;

	xt=outerop(t/hPI,xx,'*');
	fXdec=cos(xt).*repmat(rephip,1,length(xx))+sin(xt).*repmat(imphip,1,length(xx));
	fXdec=sum(fXdec.*repmat(phiK(t),1,length(xx)),1)/(2*pi)*dt/hPI;
end

function hPI = PI_deconvUestth4(W, tlim, ppphiU, hatvarU)
    phiK = @(t) (1-t.^2).^3;
    muK2 = 6;
    n = length(W);

    %grid of h values where to search for a solution
    maxh=(max(W)-min(W))/10;
    hnaive=1.06*sqrt(var(W))*n^(-1/5);
    hgrid=hnaive/3:(maxh-hnaive/3)/100:maxh;

    lh = length(hgrid);

    stdevx = max(sqrt(var(W) - hatvarU),1/n); % std(X)
    th4 = stdevx^(-9)*105/(32*sqrt(pi)); % Estimate theta4 by NR     

    dt = .0002;
    t = (-1:dt:1)';
    t=reshape(t,length(t),1);
    hgrid=reshape(hgrid,1,lh);

    toverh=t*(1./hgrid);

    phiK2=(phiK(t)).^2;
    phiU2=phiUspline(toverh,hatvarU,tlim,ppphiU).^2;

    rr=3;
    % Find h3 for th3
    term1= -hgrid.^2*muK2*th4;
    term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
    term2=sum(term2,1)*dt./(2*pi*n*hgrid.^(2*rr+1));

    ABias2 = (term1 + term2).^2;
    indh3=find(ABias2==min(ABias2),1,'first');
    h3 = hgrid(indh3);

    OO=outerop(t/h3,W,'*');
    %Estimate empirical characersitic fucntion of W
    rehatphiW=sum(cos(OO),2)/n;
    imhatphiW=sum(sin(OO),2)/n;
    clear OO;
    normhatphiW2=rehatphiW.^2+imhatphiW.^2;
    th3 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh3))*dt/(2*pi*h3^(2*rr+1));

    rr=2;
    % Find h2 for th2
    term1= -hgrid.^2*muK2*th3;
    term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
    term2=sum(term2,1)*dt./(2*pi*n*hgrid.^(2*rr+1));

    ABias2 = (term1 + term2).^2;
    indh2=find(ABias2==min(ABias2),1,'first');
    h2 = hgrid(indh2);

    OO=outerop(t/h2,W,'*');
    %Estimate empirical characersitic fucntion of W
    rehatphiW=sum(cos(OO),2)/n;
    imhatphiW=sum(sin(OO),2)/n;
    clear OO;
    normhatphiW2=rehatphiW.^2+imhatphiW.^2;
    th2 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh2))*dt/(2*pi*h2^(2*rr+1));

    term1=hgrid.^4*muK2^2*th2/4;
    term2=repmat(phiK2,1,lh)./phiU2;
    term2=sum(term2,1)*dt./(2*pi*n*hgrid);
    AMISE=term1+term2;

    indh=find(AMISE==min(AMISE),1,'first');
    hPI = hgrid(indh);
end

function y = phiUspline(t, hatvarU, tlim, ppphiU)
    ind1=(t>=tlim(1))&(t<=tlim(2));
    ind2=(t<tlim(1))|(t>tlim(2));

    phiULap = @(t) 1./(1+hatvarU/2*t.^2);


    y=0*t;
    y(ind1)=ppval(ppphiU,t(ind1));
    y(ind2)=phiULap(t(ind2));
end