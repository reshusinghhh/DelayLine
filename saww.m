length    = 14.75*1e-3;
width     = 3*1e-3;
s         = 0.25*1e-3;
thickness = 1.5*1e-3;
viaDia    = 0.5*1e-3*2;
feedWidth = viaDia/2;
gndLength = length+2*s;
gndWidth  = gndLength;
pr1 = em.internal.makerectangle(length,width)';
pr2 = em.internal.makerectangle(width,length)';

f1  = em.internal.makerectangle(feedWidth,feedWidth);

pr3 = em.internal.translateshape(f1,[-5*1e-3+feedWidth/2 0 0])';

radiator = customAntennaGeometry('Boundary',{pr1,pr2,pr3},'Operation','P1+P2+P3');

radiator.FeedLocation = [-5*1e-3 0 0];
radiator.FeedWidth = feedWidth;

ant = reflector('Exciter',radiator,'GroundPlaneLength',gndLength,'GroundPlaneWidth',gndWidth,...
    'Spacing',thickness,'EnableProbeFeed',true);

figure;
show(ant);

figure;
mesh(ant,'MaxEdgeLength',0.1);

freq = 9.14e9;
figure;impedance(ant,freq*(0.8:0.01:1.1));
