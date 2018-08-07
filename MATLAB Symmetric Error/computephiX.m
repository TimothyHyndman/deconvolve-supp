function [rephip, imphip, normphip] = computephiX(tt, theta, p)
    OO=outerop(tt,theta,'*');
    pmat=repmat(p,length(tt),1);
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