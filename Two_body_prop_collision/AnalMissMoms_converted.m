M1(1)=MX1(1)- y1 ; 
M1(2)=MX1(2)- y2 ; 
M1(3)=MX1(3)- y3 ; 
M2(1)=MX2(1)- 2*MX1(1)*y1+ y1^2 ; 
M2(2)=MX2(2)-MX1(1)*y3-MX1(3)*y1+ y1*y3 ; 
M2(3)=MX2(3)-MX1(1)*y2-MX1(2)*y1+ y1*y2 ; 
M2(4)=MX2(4)-MX1(2)*y3-MX1(3)*y2+ y2*y3 ; 
M2(5)=MX2(5)- 2*MX1(3)*y3+ y3^2 ; 
M2(6)=MX2(6)- 2*MX1(2)*y2+ y2^2 ; 
M3(1)=MX3(1)- 3*MX2(1)*y1+ 3*MX1(1)*y1^2- y1^3 ; 
M3(2)=MX3(2)+MX1(3)*y1^2-MX2(1)*y3- y1^2*y3 - 2*MX2(2)*y1+ 2*MX1(1)*y1*y3; 
M3(3)=MX3(3)+MX1(2)*y1^2-MX2(1)*y2- y1^2*y2 - 2*MX2(3)*y1+ 2*MX1(1)*y1*y2; 
M3(4)=MX3(4)-MX2(3)*y3-MX2(2)*y2-MX2(4)*y1+MX1(1)*y2*y3+MX1(2)*y1*y3+MX1(3)*y1*y2- y1*y2*y3 ; 
M3(5)=MX3(5)+MX1(1)*y3^2-MX2(5)*y1- y1*y3^2 - 2*MX2(2)*y3+ 2*MX1(3)*y1*y3; 
M3(6)=MX3(6)+MX1(1)*y2^2-MX2(6)*y1- y1*y2^2 - 2*MX2(3)*y2+ 2*MX1(2)*y1*y2; 
M3(7)=MX3(7)+MX1(3)*y2^2-MX2(6)*y3- y2^2*y3 - 2*MX2(4)*y2+ 2*MX1(2)*y2*y3; 
M3(8)=MX3(8)+MX1(2)*y3^2-MX2(5)*y2- y2*y3^2 - 2*MX2(4)*y3+ 2*MX1(3)*y2*y3; 
M3(9)=MX3(9)- 3*MX2(5)*y3+ 3*MX1(3)*y3^2- y3^3 ; 
M3(10)=MX3(10)- 3*MX2(6)*y2+ 3*MX1(2)*y2^2- y2^3 ; 
M4(1)=MX4(1)- 4*MX3(1)*y1+ 6*MX2(1)*y1^2- 4*MX1(1)*y1^3+ y1^4 ; 
M4(2)=MX4(2)-MX1(3)*y1^3-MX3(1)*y3+ y1^3*y3 + 3*MX2(2)*y1^2- 3*MX3(2)*y1- 3*MX1(1)*y1^2*y3+ 3*MX2(1)*y1*y3; 
M4(3)=MX4(3)-MX1(2)*y1^3-MX3(1)*y2+ y1^3*y2 + 3*MX2(3)*y1^2- 3*MX3(3)*y1- 3*MX1(1)*y1^2*y2+ 3*MX2(1)*y1*y2; 
M4(4)=MX4(4)+MX2(4)*y1^2-MX3(3)*y3-MX3(2)*y2-MX1(2)*y1^2*y3-MX1(3)*y1^2*y2+MX2(1)*y2*y3+ y1^2*y2*y3 - 2*MX3(4)*y1+ 2*MX2(3)*y1*y3+ 2*MX2(2)*y1*y2- 2*MX1(1)*y1*y2*y3; 
M4(5)=MX4(5)- 2*MX3(2)*y3+MX2(1)*y3^2- 2*MX3(5)*y1+ 4*MX2(2)*y1*y3- 2*MX1(1)*y1*y3^2+MX2(5)*y1^2- 2*MX1(3)*y1^2*y3+ y1^2*y3^2 ; 
M4(6)=MX4(6)- 2*MX3(3)*y2+MX2(1)*y2^2- 2*MX3(6)*y1+ 4*MX2(3)*y1*y2- 2*MX1(1)*y1*y2^2+MX2(6)*y1^2- 2*MX1(2)*y1^2*y2+ y1^2*y2^2 ; 
M4(7)=MX4(7)+MX2(2)*y2^2-MX3(6)*y3-MX3(7)*y1-MX1(1)*y2^2*y3-MX1(3)*y1*y2^2+MX2(6)*y1*y3+ y1*y2^2*y3 - 2*MX3(4)*y2+ 2*MX2(3)*y2*y3+ 2*MX2(4)*y1*y2- 2*MX1(2)*y1*y2*y3; 
M4(8)=MX4(8)+MX2(3)*y3^2-MX3(5)*y2-MX3(8)*y1-MX1(1)*y2*y3^2-MX1(2)*y1*y3^2+MX2(5)*y1*y2+ y1*y2*y3^2 - 2*MX3(4)*y3+ 2*MX2(2)*y2*y3+ 2*MX2(4)*y1*y3- 2*MX1(3)*y1*y2*y3; 
M4(9)=MX4(9)-MX1(1)*y3^3-MX3(9)*y1+ y1*y3^3 + 3*MX2(2)*y3^2- 3*MX3(5)*y3- 3*MX1(3)*y1*y3^2+ 3*MX2(5)*y1*y3; 
M4(10)=MX4(10)-MX1(1)*y2^3-MX3(10)*y1+ y1*y2^3 + 3*MX2(3)*y2^2- 3*MX3(6)*y2- 3*MX1(2)*y1*y2^2+ 3*MX2(6)*y1*y2; 
M4(11)=MX4(11)-MX1(3)*y2^3-MX3(10)*y3+ y2^3*y3 + 3*MX2(4)*y2^2- 3*MX3(7)*y2- 3*MX1(2)*y2^2*y3+ 3*MX2(6)*y2*y3; 
M4(12)=MX4(12)- 2*MX3(7)*y3+MX2(6)*y3^2- 2*MX3(8)*y2+ 4*MX2(4)*y2*y3- 2*MX1(2)*y2*y3^2+MX2(5)*y2^2- 2*MX1(3)*y2^2*y3+ y2^2*y3^2 ; 
M4(13)=MX4(13)-MX1(2)*y3^3-MX3(9)*y2+ y2*y3^3 + 3*MX2(4)*y3^2- 3*MX3(8)*y3- 3*MX1(3)*y2*y3^2+ 3*MX2(5)*y2*y3; 
M4(14)=MX4(14)- 4*MX3(9)*y3+ 6*MX2(5)*y3^2- 4*MX1(3)*y3^3+ y3^4 ; 
M4(15)=MX4(15)- 4*MX3(10)*y2+ 6*MX2(6)*y2^2- 4*MX1(2)*y2^3+ y2^4 ; 
M5(1)=MX5(1)- 5*MX4(1)*y1+ 10*MX3(1)*y1^2- 10*MX2(1)*y1^3+ 5*MX1(1)*y1^4- y1^5 ; 
M5(2)=MX5(2)+MX1(3)*y1^4-MX4(1)*y3- y1^4*y3 - 4*MX2(2)*y1^3- 4*MX4(2)*y1+ 4*MX1(1)*y1^3*y3+ 4*MX3(1)*y1*y3+ 6*MX3(2)*y1^2- 6*MX2(1)*y1^2*y3; 
M5(3)=MX5(3)+MX1(2)*y1^4-MX4(1)*y2- y1^4*y2 - 4*MX2(3)*y1^3- 4*MX4(3)*y1+ 4*MX1(1)*y1^3*y2+ 4*MX3(1)*y1*y2+ 6*MX3(3)*y1^2- 6*MX2(1)*y1^2*y2; 
M5(4)=MX5(4)-MX2(4)*y1^3-MX4(3)*y3-MX4(2)*y2+MX1(2)*y1^3*y3+MX1(3)*y1^3*y2+MX3(1)*y2*y3- y1^3*y2*y3 + 3*MX3(4)*y1^2- 3*MX4(4)*y1- 3*MX2(3)*y1^2*y3- 3*MX2(2)*y1^2*y2+ 3*MX3(3)*y1*y3+ 3*MX3(2)*y1*y2+ 3*MX1(1)*y1^2*y2*y3- 3*MX2(1)*y1*y2*y3; 
M5(5)=MX5(5)- 2*MX4(2)*y3+MX3(1)*y3^2- 3*MX4(5)*y1+ 6*MX3(2)*y1*y3- 3*MX2(1)*y1*y3^2+ 3*MX3(5)*y1^2- 6*MX2(2)*y1^2*y3+ 3*MX1(1)*y1^2*y3^2-MX2(5)*y1^3+ 2*MX1(3)*y1^3*y3- y1^3*y3^2 ; 
M5(6)=MX5(6)- 2*MX4(3)*y2+MX3(1)*y2^2- 3*MX4(6)*y1+ 6*MX3(3)*y1*y2- 3*MX2(1)*y1*y2^2+ 3*MX3(6)*y1^2- 6*MX2(3)*y1^2*y2+ 3*MX1(1)*y1^2*y2^2-MX2(6)*y1^3+ 2*MX1(2)*y1^3*y2- y1^3*y2^2 ; 
M5(7)=MX5(7)+MX3(2)*y2^2-MX4(6)*y3+MX3(7)*y1^2+MX1(3)*y1^2*y2^2-MX2(1)*y2^2*y3-MX2(6)*y1^2*y3- y1^2*y2^2*y3 - 2*MX4(7)*y1- 2*MX4(4)*y2- 2*MX2(2)*y1*y2^2+ 2*MX3(6)*y1*y3- 2*MX2(4)*y1^2*y2+ 2*MX3(3)*y2*y3+ 2*MX1(1)*y1*y2^2*y3+ 2*MX1(2)*y1^2*y2*y3+ 4*MX3(4)*y1*y2- 4*MX2(3)*y1*y2*y3; 
M5(8)=MX5(8)+MX3(8)*y1^2+MX3(3)*y3^2-MX4(5)*y2+MX1(2)*y1^2*y3^2-MX2(1)*y2*y3^2-MX2(5)*y1^2*y2- y1^2*y2*y3^2 - 2*MX4(8)*y1- 2*MX4(4)*y3- 2*MX2(3)*y1*y3^2+ 2*MX3(5)*y1*y2- 2*MX2(4)*y1^2*y3+ 2*MX3(2)*y2*y3+ 2*MX1(1)*y1*y2*y3^2+ 2*MX1(3)*y1^2*y2*y3+ 4*MX3(4)*y1*y3- 4*MX2(2)*y1*y2*y3; 
M5(9)=MX5(9)- 3*MX4(5)*y3+ 3*MX3(2)*y3^2-MX2(1)*y3^3- 2*MX4(9)*y1+ 6*MX3(5)*y1*y3- 6*MX2(2)*y1*y3^2+ 2*MX1(1)*y1*y3^3+MX3(9)*y1^2- 3*MX2(5)*y1^2*y3+ 3*MX1(3)*y1^2*y3^2- y1^2*y3^3 ; 
M5(10)=MX5(10)- 3*MX4(6)*y2+ 3*MX3(3)*y2^2-MX2(1)*y2^3- 2*MX4(10)*y1+ 6*MX3(6)*y1*y2- 6*MX2(3)*y1*y2^2+ 2*MX1(1)*y1*y2^3+MX3(10)*y1^2- 3*MX2(6)*y1^2*y2+ 3*MX1(2)*y1^2*y2^2- y1^2*y2^3 ; 
M5(11)=MX5(11)-MX2(2)*y2^3-MX4(10)*y3-MX4(11)*y1+MX1(1)*y2^3*y3+MX1(3)*y1*y2^3+MX3(10)*y1*y3- y1*y2^3*y3 + 3*MX3(4)*y2^2- 3*MX4(7)*y2- 3*MX2(3)*y2^2*y3+ 3*MX3(6)*y2*y3- 3*MX2(4)*y1*y2^2+ 3*MX3(7)*y1*y2+ 3*MX1(2)*y1*y2^2*y3- 3*MX2(6)*y1*y2*y3; 
M5(12)=MX5(12)+MX3(6)*y3^2+MX3(5)*y2^2-MX4(12)*y1+MX1(1)*y2^2*y3^2-MX2(6)*y1*y3^2-MX2(5)*y1*y2^2- y1*y2^2*y3^2 - 2*MX4(8)*y2- 2*MX4(7)*y3- 2*MX2(3)*y2*y3^2+ 2*MX3(8)*y1*y2- 2*MX2(2)*y2^2*y3+ 2*MX3(7)*y1*y3+ 2*MX1(2)*y1*y2*y3^2+ 2*MX1(3)*y1*y2^2*y3+ 4*MX3(4)*y2*y3- 4*MX2(4)*y1*y2*y3; 
M5(13)=MX5(13)-MX2(3)*y3^3-MX4(9)*y2-MX4(13)*y1+MX1(1)*y2*y3^3+MX1(2)*y1*y3^3+MX3(9)*y1*y2- y1*y2*y3^3 + 3*MX3(4)*y3^2- 3*MX4(8)*y3- 3*MX2(2)*y2*y3^2+ 3*MX3(5)*y2*y3- 3*MX2(4)*y1*y3^2+ 3*MX3(8)*y1*y3+ 3*MX1(3)*y1*y2*y3^2- 3*MX2(5)*y1*y2*y3; 
M5(14)=MX5(14)+MX1(1)*y3^4-MX4(14)*y1- y1*y3^4 - 4*MX2(2)*y3^3- 4*MX4(9)*y3+ 4*MX1(3)*y1*y3^3+ 4*MX3(9)*y1*y3+ 6*MX3(5)*y3^2- 6*MX2(5)*y1*y3^2; 
M5(15)=MX5(15)+MX1(1)*y2^4-MX4(15)*y1- y1*y2^4 - 4*MX2(3)*y2^3- 4*MX4(10)*y2+ 4*MX1(2)*y1*y2^3+ 4*MX3(10)*y1*y2+ 6*MX3(6)*y2^2- 6*MX2(6)*y1*y2^2; 
M5(16)=MX5(16)+MX1(3)*y2^4-MX4(15)*y3- y2^4*y3 - 4*MX2(4)*y2^3- 4*MX4(11)*y2+ 4*MX1(2)*y2^3*y3+ 4*MX3(10)*y2*y3+ 6*MX3(7)*y2^2- 6*MX2(6)*y2^2*y3; 
M5(17)=MX5(17)- 2*MX4(11)*y3+MX3(10)*y3^2- 3*MX4(12)*y2+ 6*MX3(7)*y2*y3- 3*MX2(6)*y2*y3^2+ 3*MX3(8)*y2^2- 6*MX2(4)*y2^2*y3+ 3*MX1(2)*y2^2*y3^2-MX2(5)*y2^3+ 2*MX1(3)*y2^3*y3- y2^3*y3^2 ; 
M5(18)=MX5(18)- 3*MX4(12)*y3+ 3*MX3(7)*y3^2-MX2(6)*y3^3- 2*MX4(13)*y2+ 6*MX3(8)*y2*y3- 6*MX2(4)*y2*y3^2+ 2*MX1(2)*y2*y3^3+MX3(9)*y2^2- 3*MX2(5)*y2^2*y3+ 3*MX1(3)*y2^2*y3^2- y2^2*y3^3 ; 
M5(19)=MX5(19)+MX1(2)*y3^4-MX4(14)*y2- y2*y3^4 - 4*MX2(4)*y3^3- 4*MX4(13)*y3+ 4*MX1(3)*y2*y3^3+ 4*MX3(9)*y2*y3+ 6*MX3(8)*y3^2- 6*MX2(5)*y2*y3^2; 
M5(20)=MX5(20)- 5*MX4(14)*y3+ 10*MX3(9)*y3^2- 10*MX2(5)*y3^3+ 5*MX1(3)*y3^4- y3^5 ; 
M5(21)=MX5(21)- 5*MX4(15)*y2+ 10*MX3(10)*y2^2- 10*MX2(6)*y2^3+ 5*MX1(2)*y2^4- y2^5 ; 
M6(1)=MX6(1)- 6*MX5(1)*y1+ 15*MX4(1)*y1^2- 20*MX3(1)*y1^3+ 15*MX2(1)*y1^4- 6*MX1(1)*y1^5+ y1^6 ; 
M6(2)=MX6(2)-MX1(3)*y1^5-MX5(1)*y3+ y1^5*y3 + 5*MX2(2)*y1^4- 5*MX5(2)*y1- 5*MX1(1)*y1^4*y3+ 5*MX4(1)*y1*y3- 10*MX3(2)*y1^3+ 10*MX4(2)*y1^2+ 10*MX2(1)*y1^3*y3- 10*MX3(1)*y1^2*y3; 
M6(3)=MX6(3)-MX1(2)*y1^5-MX5(1)*y2+ y1^5*y2 + 5*MX2(3)*y1^4- 5*MX5(3)*y1- 5*MX1(1)*y1^4*y2+ 5*MX4(1)*y1*y2- 10*MX3(3)*y1^3+ 10*MX4(3)*y1^2+ 10*MX2(1)*y1^3*y2- 10*MX3(1)*y1^2*y2; 
M6(4)=MX6(4)+MX2(4)*y1^4-MX5(3)*y3-MX5(2)*y2-MX1(2)*y1^4*y3-MX1(3)*y1^4*y2+MX4(1)*y2*y3+ y1^4*y2*y3 + 6*MX4(4)*y1^2- 6*MX3(3)*y1^2*y3- 6*MX3(2)*y1^2*y2+ 6*MX2(1)*y1^2*y2*y3- 4*MX3(4)*y1^3- 4*MX5(4)*y1+ 4*MX2(3)*y1^3*y3+ 4*MX2(2)*y1^3*y2+ 4*MX4(3)*y1*y3+ 4*MX4(2)*y1*y2- 4*MX1(1)*y1^3*y2*y3- 4*MX3(1)*y1*y2*y3; 
M6(5)=MX6(5)- 2*MX5(2)*y3+MX4(1)*y3^2- 4*MX5(5)*y1+ 8*MX4(2)*y1*y3- 4*MX3(1)*y1*y3^2+ 6*MX4(5)*y1^2- 12*MX3(2)*y1^2*y3+ 6*MX2(1)*y1^2*y3^2- 4*MX3(5)*y1^3+ 8*MX2(2)*y1^3*y3- 4*MX1(1)*y1^3*y3^2+MX2(5)*y1^4- 2*MX1(3)*y1^4*y3+ y1^4*y3^2 ; 
M6(6)=MX6(6)- 2*MX5(3)*y2+MX4(1)*y2^2- 4*MX5(6)*y1+ 8*MX4(3)*y1*y2- 4*MX3(1)*y1*y2^2+ 6*MX4(6)*y1^2- 12*MX3(3)*y1^2*y2+ 6*MX2(1)*y1^2*y2^2- 4*MX3(6)*y1^3+ 8*MX2(3)*y1^3*y2- 4*MX1(1)*y1^3*y2^2+MX2(6)*y1^4- 2*MX1(2)*y1^4*y2+ y1^4*y2^2 ; 
M6(7)=MX6(7)+MX4(2)*y2^2-MX5(6)*y3-MX3(7)*y1^3-MX1(3)*y1^3*y2^2-MX3(1)*y2^2*y3+MX2(6)*y1^3*y3+ y1^3*y2^2*y3 + 3*MX4(7)*y1^2- 3*MX5(7)*y1+ 3*MX2(2)*y1^2*y2^2- 3*MX3(6)*y1^2*y3- 3*MX3(2)*y1*y2^2+ 3*MX4(6)*y1*y3- 3*MX1(1)*y1^2*y2^2*y3+ 3*MX2(1)*y1*y2^2*y3- 2*MX5(4)*y2+ 2*MX2(4)*y1^3*y2+ 2*MX4(3)*y2*y3- 2*MX1(2)*y1^3*y2*y3- 6*MX3(4)*y1^2*y2+ 6*MX4(4)*y1*y2+ 6*MX2(3)*y1^2*y2*y3- 6*MX3(3)*y1*y2*y3; 
M6(8)=MX6(8)-MX3(8)*y1^3+MX4(3)*y3^2-MX5(5)*y2-MX1(2)*y1^3*y3^2-MX3(1)*y2*y3^2+MX2(5)*y1^3*y2+ y1^3*y2*y3^2 + 3*MX4(8)*y1^2- 3*MX5(8)*y1+ 3*MX2(3)*y1^2*y3^2- 3*MX3(5)*y1^2*y2- 3*MX3(3)*y1*y3^2+ 3*MX4(5)*y1*y2- 3*MX1(1)*y1^2*y2*y3^2+ 3*MX2(1)*y1*y2*y3^2- 2*MX5(4)*y3+ 2*MX2(4)*y1^3*y3+ 2*MX4(2)*y2*y3- 2*MX1(3)*y1^3*y2*y3- 6*MX3(4)*y1^2*y3+ 6*MX4(4)*y1*y3+ 6*MX2(2)*y1^2*y2*y3- 6*MX3(2)*y1*y2*y3; 
M6(9)=MX6(9)- 3*MX5(5)*y3+ 3*MX4(2)*y3^2-MX3(1)*y3^3- 3*MX5(9)*y1+ 9*MX4(5)*y1*y3- 9*MX3(2)*y1*y3^2+ 3*MX2(1)*y1*y3^3+ 3*MX4(9)*y1^2- 9*MX3(5)*y1^2*y3+ 9*MX2(2)*y1^2*y3^2- 3*MX1(1)*y1^2*y3^3-MX3(9)*y1^3+ 3*MX2(5)*y1^3*y3- 3*MX1(3)*y1^3*y3^2+ y1^3*y3^3 ; 
M6(10)=MX6(10)- 3*MX5(6)*y2+ 3*MX4(3)*y2^2-MX3(1)*y2^3- 3*MX5(10)*y1+ 9*MX4(6)*y1*y2- 9*MX3(3)*y1*y2^2+ 3*MX2(1)*y1*y2^3+ 3*MX4(10)*y1^2- 9*MX3(6)*y1^2*y2+ 9*MX2(3)*y1^2*y2^2- 3*MX1(1)*y1^2*y2^3-MX3(10)*y1^3+ 3*MX2(6)*y1^3*y2- 3*MX1(2)*y1^3*y2^2+ y1^3*y2^3 ; 
M6(11)=MX6(11)-MX3(2)*y2^3-MX5(10)*y3+MX4(11)*y1^2-MX1(3)*y1^2*y2^3+MX2(1)*y2^3*y3-MX3(10)*y1^2*y3+ y1^2*y2^3*y3 + 3*MX4(4)*y2^2- 3*MX5(7)*y2+ 3*MX2(4)*y1^2*y2^2- 3*MX3(3)*y2^2*y3+ 3*MX4(6)*y2*y3- 3*MX3(7)*y1^2*y2- 3*MX1(2)*y1^2*y2^2*y3+ 3*MX2(6)*y1^2*y2*y3- 2*MX5(11)*y1+ 2*MX2(2)*y1*y2^3+ 2*MX4(10)*y1*y3- 2*MX1(1)*y1*y2^3*y3- 6*MX3(4)*y1*y2^2+ 6*MX4(7)*y1*y2+ 6*MX2(3)*y1*y2^2*y3- 6*MX3(6)*y1*y2*y3; 
M6(12)=MX6(12)- 2*MX5(7)*y3+MX4(6)*y3^2- 2*MX5(8)*y2+ 4*MX4(4)*y2*y3- 2*MX3(3)*y2*y3^2+MX4(5)*y2^2- 2*MX3(2)*y2^2*y3+MX2(1)*y2^2*y3^2- 2*MX5(12)*y1+ 4*MX4(7)*y1*y3- 2*MX3(6)*y1*y3^2+ 4*MX4(8)*y1*y2- 8*MX3(4)*y1*y2*y3+ 4*MX2(3)*y1*y2*y3^2- 2*MX3(5)*y1*y2^2+ 4*MX2(2)*y1*y2^2*y3- 2*MX1(1)*y1*y2^2*y3^2+MX4(12)*y1^2- 2*MX3(7)*y1^2*y3+MX2(6)*y1^2*y3^2- 2*MX3(8)*y1^2*y2+ 4*MX2(4)*y1^2*y2*y3- 2*MX1(2)*y1^2*y2*y3^2+MX2(5)*y1^2*y2^2- 2*MX1(3)*y1^2*y2^2*y3+ y1^2*y2^2*y3^2 ; 
M6(13)=MX6(13)+MX4(13)*y1^2-MX3(3)*y3^3-MX5(9)*y2-MX1(2)*y1^2*y3^3+MX2(1)*y2*y3^3-MX3(9)*y1^2*y2+ y1^2*y2*y3^3 + 3*MX4(4)*y3^2- 3*MX5(8)*y3+ 3*MX2(4)*y1^2*y3^2- 3*MX3(8)*y1^2*y3- 3*MX3(2)*y2*y3^2+ 3*MX4(5)*y2*y3- 3*MX1(3)*y1^2*y2*y3^2+ 3*MX2(5)*y1^2*y2*y3- 2*MX5(13)*y1+ 2*MX2(3)*y1*y3^3+ 2*MX4(9)*y1*y2- 2*MX1(1)*y1*y2*y3^3- 6*MX3(4)*y1*y3^2+ 6*MX4(8)*y1*y3+ 6*MX2(2)*y1*y2*y3^2- 6*MX3(5)*y1*y2*y3; 
M6(14)=MX6(14)- 4*MX5(9)*y3+ 6*MX4(5)*y3^2- 4*MX3(2)*y3^3+MX2(1)*y3^4- 2*MX5(14)*y1+ 8*MX4(9)*y1*y3- 12*MX3(5)*y1*y3^2+ 8*MX2(2)*y1*y3^3- 2*MX1(1)*y1*y3^4+MX4(14)*y1^2- 4*MX3(9)*y1^2*y3+ 6*MX2(5)*y1^2*y3^2- 4*MX1(3)*y1^2*y3^3+ y1^2*y3^4 ; 
M6(15)=MX6(15)- 4*MX5(10)*y2+ 6*MX4(6)*y2^2- 4*MX3(3)*y2^3+MX2(1)*y2^4- 2*MX5(15)*y1+ 8*MX4(10)*y1*y2- 12*MX3(6)*y1*y2^2+ 8*MX2(3)*y1*y2^3- 2*MX1(1)*y1*y2^4+MX4(15)*y1^2- 4*MX3(10)*y1^2*y2+ 6*MX2(6)*y1^2*y2^2- 4*MX1(2)*y1^2*y2^3+ y1^2*y2^4 ; 
M6(16)=MX6(16)+MX2(2)*y2^4-MX5(15)*y3-MX5(16)*y1-MX1(1)*y2^4*y3-MX1(3)*y1*y2^4+MX4(15)*y1*y3+ y1*y2^4*y3 + 6*MX4(7)*y2^2- 6*MX3(6)*y2^2*y3- 6*MX3(7)*y1*y2^2+ 6*MX2(6)*y1*y2^2*y3- 4*MX3(4)*y2^3- 4*MX5(11)*y2+ 4*MX2(3)*y2^3*y3+ 4*MX4(10)*y2*y3+ 4*MX2(4)*y1*y2^3+ 4*MX4(11)*y1*y2- 4*MX1(2)*y1*y2^3*y3- 4*MX3(10)*y1*y2*y3; 
M6(17)=MX6(17)+MX4(10)*y3^2-MX3(5)*y2^3-MX5(17)*y1-MX1(1)*y2^3*y3^2-MX3(10)*y1*y3^2+MX2(5)*y1*y2^3+ y1*y2^3*y3^2 + 3*MX4(8)*y2^2- 3*MX5(12)*y2+ 3*MX2(3)*y2^2*y3^2- 3*MX3(6)*y2*y3^2- 3*MX3(8)*y1*y2^2+ 3*MX4(12)*y1*y2- 3*MX1(2)*y1*y2^2*y3^2+ 3*MX2(6)*y1*y2*y3^2- 2*MX5(11)*y3+ 2*MX2(2)*y2^3*y3+ 2*MX4(11)*y1*y3- 2*MX1(3)*y1*y2^3*y3- 6*MX3(4)*y2^2*y3+ 6*MX4(7)*y2*y3+ 6*MX2(4)*y1*y2^2*y3- 6*MX3(7)*y1*y2*y3; 
M6(18)=MX6(18)-MX3(6)*y3^3+MX4(9)*y2^2-MX5(18)*y1-MX1(1)*y2^2*y3^3+MX2(6)*y1*y3^3-MX3(9)*y1*y2^2+ y1*y2^2*y3^3 + 3*MX4(7)*y3^2- 3*MX5(12)*y3+ 3*MX2(2)*y2^2*y3^2- 3*MX3(5)*y2^2*y3- 3*MX3(7)*y1*y3^2+ 3*MX4(12)*y1*y3- 3*MX1(3)*y1*y2^2*y3^2+ 3*MX2(5)*y1*y2^2*y3- 2*MX5(13)*y2+ 2*MX2(3)*y2*y3^3+ 2*MX4(13)*y1*y2- 2*MX1(2)*y1*y2*y3^3- 6*MX3(4)*y2*y3^2+ 6*MX4(8)*y2*y3+ 6*MX2(4)*y1*y2*y3^2- 6*MX3(8)*y1*y2*y3; 
M6(19)=MX6(19)+MX2(3)*y3^4-MX5(14)*y2-MX5(19)*y1-MX1(1)*y2*y3^4-MX1(2)*y1*y3^4+MX4(14)*y1*y2+ y1*y2*y3^4 + 6*MX4(8)*y3^2- 6*MX3(5)*y2*y3^2- 6*MX3(8)*y1*y3^2+ 6*MX2(5)*y1*y2*y3^2- 4*MX3(4)*y3^3- 4*MX5(13)*y3+ 4*MX2(2)*y2*y3^3+ 4*MX4(9)*y2*y3+ 4*MX2(4)*y1*y3^3+ 4*MX4(13)*y1*y3- 4*MX1(3)*y1*y2*y3^3- 4*MX3(9)*y1*y2*y3; 
M6(20)=MX6(20)-MX1(1)*y3^5-MX5(20)*y1+ y1*y3^5 + 5*MX2(2)*y3^4- 5*MX5(14)*y3- 5*MX1(3)*y1*y3^4+ 5*MX4(14)*y1*y3- 10*MX3(5)*y3^3+ 10*MX4(9)*y3^2+ 10*MX2(5)*y1*y3^3- 10*MX3(9)*y1*y3^2; 
M6(21)=MX6(21)-MX1(1)*y2^5-MX5(21)*y1+ y1*y2^5 + 5*MX2(3)*y2^4- 5*MX5(15)*y2- 5*MX1(2)*y1*y2^4+ 5*MX4(15)*y1*y2- 10*MX3(6)*y2^3+ 10*MX4(10)*y2^2+ 10*MX2(6)*y1*y2^3- 10*MX3(10)*y1*y2^2; 
M6(22)=MX6(22)-MX1(3)*y2^5-MX5(21)*y3+ y2^5*y3 + 5*MX2(4)*y2^4- 5*MX5(16)*y2- 5*MX1(2)*y2^4*y3+ 5*MX4(15)*y2*y3- 10*MX3(7)*y2^3+ 10*MX4(11)*y2^2+ 10*MX2(6)*y2^3*y3- 10*MX3(10)*y2^2*y3; 
M6(23)=MX6(23)- 2*MX5(16)*y3+MX4(15)*y3^2- 4*MX5(17)*y2+ 8*MX4(11)*y2*y3- 4*MX3(10)*y2*y3^2+ 6*MX4(12)*y2^2- 12*MX3(7)*y2^2*y3+ 6*MX2(6)*y2^2*y3^2- 4*MX3(8)*y2^3+ 8*MX2(4)*y2^3*y3- 4*MX1(2)*y2^3*y3^2+MX2(5)*y2^4- 2*MX1(3)*y2^4*y3+ y2^4*y3^2 ; 
M6(24)=MX6(24)- 3*MX5(17)*y3+ 3*MX4(11)*y3^2-MX3(10)*y3^3- 3*MX5(18)*y2+ 9*MX4(12)*y2*y3- 9*MX3(7)*y2*y3^2+ 3*MX2(6)*y2*y3^3+ 3*MX4(13)*y2^2- 9*MX3(8)*y2^2*y3+ 9*MX2(4)*y2^2*y3^2- 3*MX1(2)*y2^2*y3^3-MX3(9)*y2^3+ 3*MX2(5)*y2^3*y3- 3*MX1(3)*y2^3*y3^2+ y2^3*y3^3 ; 
M6(25)=MX6(25)- 4*MX5(18)*y3+ 6*MX4(12)*y3^2- 4*MX3(7)*y3^3+MX2(6)*y3^4- 2*MX5(19)*y2+ 8*MX4(13)*y2*y3- 12*MX3(8)*y2*y3^2+ 8*MX2(4)*y2*y3^3- 2*MX1(2)*y2*y3^4+MX4(14)*y2^2- 4*MX3(9)*y2^2*y3+ 6*MX2(5)*y2^2*y3^2- 4*MX1(3)*y2^2*y3^3+ y2^2*y3^4 ; 
M6(26)=MX6(26)-MX1(2)*y3^5-MX5(20)*y2+ y2*y3^5 + 5*MX2(4)*y3^4- 5*MX5(19)*y3- 5*MX1(3)*y2*y3^4+ 5*MX4(14)*y2*y3- 10*MX3(8)*y3^3+ 10*MX4(13)*y3^2+ 10*MX2(5)*y2*y3^3- 10*MX3(9)*y2*y3^2; 
M6(27)=MX6(27)- 6*MX5(20)*y3+ 15*MX4(14)*y3^2- 20*MX3(9)*y3^3+ 15*MX2(5)*y3^4- 6*MX1(3)*y3^5+ y3^6 ; 
M6(28)=MX6(28)- 6*MX5(21)*y2+ 15*MX4(15)*y2^2- 20*MX3(10)*y2^3+ 15*MX2(6)*y2^4- 6*MX1(2)*y2^5+ y2^6 ; 
 
 
 
 
 
 
 
