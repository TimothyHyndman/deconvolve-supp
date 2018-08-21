function [rehatpsi,imhatpsi,sqabshatpsi]=compute_psi_W(tt,W)

n = length(W);

WW = outerop(W,W,'-');
WW = WW(~logical(eye(n)));

rehatpsi = 0*tt;
imhatpsi = 0*tt;

for i=1:length(tt)
	rehatpsi(i) = sum(cos(tt(i)*WW));
	imhatpsi(i) = sum(sin(tt(i)*WW));
end

rehatpsi = rehatpsi' / (n * (n - 1));
imhatpsi = imhatpsi' / (n * (n - 1));

sqabshatpsi = sqrt(sqrt(rehatpsi.^2 + imhatpsi.^2));