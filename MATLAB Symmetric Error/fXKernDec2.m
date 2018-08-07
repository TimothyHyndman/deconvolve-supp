function fXdec = fXKernDec2(xx, hPI, W, tlim, ppphiU, var_U)
	phiK = @(t) (1-t.^2).^3;
	muK2 = 6;

	dt = .0002;
	t = (-1:dt:1)';
	t=reshape(t,length(t),1);
	
	phiU=phiUspline(t/hPI,var_U,tlim,ppphiU);
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

function y = phiUspline(t, var_U, tlim, ppphiU)
    ind1=(t>=tlim(1))&(t<=tlim(2));
    ind2=(t<tlim(1))|(t>tlim(2));

    phiULap = @(t) 1./(1+var_U/2*t.^2);


    y=0*t;
    y(ind1)=ppval(ppphiU,t(ind1));
    y(ind2)=phiULap(t(ind2));
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