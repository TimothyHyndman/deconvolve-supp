function hPI = PI_deconvUestth4(W, tlim, ppphiU, var_U)
    phiK = @(t) (1-t.^2).^3;
    muK2 = 6;
    RK = 1024 / 3003 / pi;
    n = length(W);

    %grid of h values where to search for a solution
    maxh=(max(W)-min(W))/10;
    hnaive = ((8 * sqrt(pi) * RK/3/muK2^2)^0.2) * sqrt(var(W))*n^(-1/5);
    hgrid=hnaive/3:(maxh-hnaive/3)/100:maxh;

    lh = length(hgrid);

    stdevx = max(sqrt(var(W) - var_U), 1 / n); % std(X)
    th4 = stdevx^(-9)*105/(32*sqrt(pi)); % Estimate theta4 by NR     

    dt = .0002;
    t = (-1:dt:1)';
    t=reshape(t,length(t),1);
    hgrid=reshape(hgrid,1,lh);

    toverh=t*(1./hgrid);

    phiK2=(phiK(t)).^2;
    phiU2=phiUspline(toverh,var_U,tlim,ppphiU).^2;

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

function y = phiUspline(t, var_U, tlim, ppphiU)
    ind1=(t>=tlim(1))&(t<=tlim(2));
    ind2=(t<tlim(1))|(t>tlim(2));

    phiULap = @(t) 1./(1+var_U/2*t.^2);


    y=0*t;
    y(ind1)=ppval(ppphiU,t(ind1));
    y(ind2)=phiULap(t(ind2));
end