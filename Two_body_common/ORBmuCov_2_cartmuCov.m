function [mux,Px]=ORBmuCov_2_cartmuCov(muorb,Porb)
[x,w]=GH_points(muorb(:),Porb,4);
XY = OE2XYZ_multiple(x);
[mux,Px]=MeanCov(XY,w);
end