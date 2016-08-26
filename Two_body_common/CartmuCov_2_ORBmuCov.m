function [muorb,Porb]=CartmuCov_2_ORBmuCov(mux,Px)
[x,w]=GH_points(mux(:),Px,4);
OO = XYZ2OE_multiple(x);
[muorb,Porb]=MeanCov(OO,w);
end