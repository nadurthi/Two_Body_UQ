load('OBEcut6')
M1cut6=M1;
M2cut6=M2;
M3cut6=M3;
M4cut6=M4;

load('OBEcut8')
M1cut8=M1;
M2cut8=M2;
M3cut8=M3;
M4cut8=M4;

load('OBEgh4')
M1gh4=M1;
M2gh4=M2;
M3gh4=M3;
M4gh4=M4;


load('OBEgh4')
M1gh5=M1;
M2gh5=M2;
M3gh5=M3;
M4gh5=M4;

load('OBEgh4')
M1gh6=M1;
M2gh6=M2;
M3gh6=M3;
M4gh6=M4;

load('OBEgh7')
M1gh7=M1;
M2gh7=M2;
M3gh7=M3;
M4gh7=M4;

load('OBEMCmoms5x100000')
M1mc=M1;
M2mc=M2;
M3mc=M3;
M4mc=M4;



disp([' cut6 ',' cut8 ',' gh4 ',' gh5 ',' gh6 ',' mc '])
[sqrt(sum(sqrt(sum((M1cut6-M1gh7).^2,2)/6).^2)/110),sqrt(sum(sqrt(sum((M1cut8-M1gh7).^2,2)/6).^2)/110),sqrt(sum(sqrt(sum((M1gh4-M1gh7).^2,2)/6).^2)/110),sqrt(sum(sqrt(sum((M1gh5-M1gh7).^2,2)/6).^2)/110),sqrt(sum(sqrt(sum((M1gh6-M1gh7).^2,2)/6).^2)/110),sqrt(sum(sqrt(sum((M1mc-M1gh7).^2,2)/6).^2)/110)]
sqrt([sqrt(sum(sqrt(sum((M2cut6-M2gh7).^2,2)/21).^2)/110),sqrt(sum(sqrt(sum((M2cut8-M2gh7).^2,2)/21).^2)/110),sqrt(sum(sqrt(sum((M2gh4-M2gh7).^2,2)/21).^2)/110),sqrt(sum(sqrt(sum((M2gh5-M2gh7).^2,2)/21).^2)/110),sqrt(sum(sqrt(sum((M2gh6-M2gh7).^2,2)/21).^2)/110),sqrt(sum(sqrt(sum((M2mc-M2gh7).^2,2)/21).^2)/110)])
([sqrt(sum(sqrt(sum((M3cut6-M3gh7).^2,2)/56).^2)/110),sqrt(sum(sqrt(sum((M3cut8-M3gh7).^2,2)/56).^2)/110),sqrt(sum(sqrt(sum((M3gh4-M3gh7).^2,2)/56).^2)/110),sqrt(sum(sqrt(sum((M3gh5-M3gh7).^2,2)/56).^2)/110),sqrt(sum(sqrt(sum((M3gh6-M3gh7).^2,2)/56).^2)/110),sqrt(sum(sqrt(sum((M3mc-M3gh7).^2,2)/56).^2)/110)]).^(1/3)
([sqrt(sum(sqrt(sum((M4cut6-M4gh7).^2,2)/126).^2)/110),sqrt(sum(sqrt(sum((M4cut8-M4gh7).^2,2)/126).^2)/110),sqrt(sum(sqrt(sum((M4gh4-M4gh7).^2,2)/126).^2)/110),sqrt(sum(sqrt(sum((M4gh5-M4gh7).^2,2)/126).^2)/110),sqrt(sum(sqrt(sum((M4gh6-M4gh7).^2,2)/126).^2)/110),sqrt(sum(sqrt(sum((M4mc-M4gh7).^2,2)/126).^2)/110)]).^(1/4)

M2cut4_equinoc
disp([' ut ',' cut4 ',' cut6 ',' cut8 ',' gh4 ',' gh5 ',' gh6 ',' mc '])
[matrixRMSEnorm(M1ut_equinoc-M1gh7_equinoc),matrixRMSEnorm(M1cut4_equinoc-M1gh7_equinoc),matrixRMSEnorm(M1cut6_equinoc-M1gh7_equinoc),matrixRMSEnorm(M1cut8_equinoc-M1gh7_equinoc),matrixRMSEnorm(M1gh4_equinoc-M1gh7_equinoc),matrixRMSEnorm(M1gh5_equinoc-M1gh7_equinoc),matrixRMSEnorm(M1gh6_equinoc-M1gh7_equinoc)]

disp([' ut ',' cut4 ',' cut6 ',' cut8 ',' gh4 ',' gh5 ',' gh6 ',' mc '])
sqrt([matrixRMSEnorm(M2ut_equinoc-M2gh7_equinoc),matrixRMSEnorm(M2cut4_equinoc-M2gh7_equinoc),matrixRMSEnorm(M2cut6_equinoc-M2gh7_equinoc),matrixRMSEnorm(M2cut8_equinoc-M2gh7_equinoc),matrixRMSEnorm(M2gh4_equinoc-M2gh7_equinoc),matrixRMSEnorm(M2gh5_equinoc-M2gh7_equinoc),matrixRMSEnorm(M2gh6_equinoc-M2gh7_equinoc)])
