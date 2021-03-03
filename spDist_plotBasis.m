% spDist_plotBasis.m
%
% plots the basis set for Fig. 2A

angs = 0:360; % in polar ang

nchans = 8;

chan_centers = linspace(360/nchans,360,nchans);


basis_set = spDist_makeX1(angs,chan_centers);

basis_colors = hsv(nchans);

%% normal plot

figure('name','Figure1B');
hold on;
for cc = 1:nchans
    plot(angs,basis_set(:,cc),'-','Color',basis_colors(cc,:));
end

xlim([angs(1) angs(end)]);
ylim([0 1]);
set(gca,'XTick',0:90:360,'YTick',[0 0.5 1],'TickDir','out');
xlabel('Polar angle (\circ)');
title('Basis set');

hold off;

