%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------



function y=phiUspline(t,hatvarU,tlim,ppphiU)

ind1=(t>=tlim(1))&(t<=tlim(2));
ind2=(t<tlim(1))|(t>tlim(2));

phiULap = @(t) 1./(1+hatvarU/2*t.^2);


y=0*t;
y(ind1)=ppval(ppphiU,t(ind1));
y(ind2)=phiULap(t(ind2));
y;
