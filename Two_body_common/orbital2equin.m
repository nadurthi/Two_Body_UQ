function [a,P1,P2,Q1,Q2,l]=orbital2equin(a,e,i,omg,Omg,M)
wb=Omg+omg;
P1=e*sin(wb);
P2=e*cos(wb);
Q1=tan(i/2)*sin(Omg);
Q2=tan(i/2)*cos(Omg);
a=a;
l=wb+M;



