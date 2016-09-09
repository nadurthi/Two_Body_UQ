function [t,y]=twoBodyKeplerProp(tvec,x)

r0=x(1:3);
v0=x(4:6);
y=zeros(length(tvec),6);
for i=1:1:length(tvec)
[r,v] = keplerUniversal(r0(:),v0(:),tvec(i)-tvec(1),398601.2);
y(i,:)=[r(:)',v(:)'];
end
t=tvec(:);

