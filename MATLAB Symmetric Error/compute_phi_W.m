function  [rehatphiW, imhatphiW] = compute_phi_W(tt, W)
    % Estimate empirical characersitic fucntion of W----------------------------
    n = length(W);
    OO = reshape(tt, [], 1) * reshape(W, 1, []);
    % OO=outerop(tt,W,'*');
    rehatphiW=sum(cos(OO),2)/n;
    imhatphiW=sum(sin(OO),2)/n;
end