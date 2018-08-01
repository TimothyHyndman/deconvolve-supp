function plot_decon_err_sym(truedens)
	%Plot stuff
	[fW,xW]=ksdensity(W,xx);

	
	hX = plot(xW,truedens(xW),'g','LineWidth',2);
	hold on
	hW = plot(xW,fW,'b','LineWidth',2);
	hj = scatter(xj,pj,'m','filled');
	plot(xW,fXdeconvoluted,'m','LineWidth',2)
	hold off
	legend([hX hW hj],{'f_X','f_W','Our estimate for f_X'},'FontSize',20,'FontWeight','bold')
	axis([min(xW) max(xW) 0 1])
end