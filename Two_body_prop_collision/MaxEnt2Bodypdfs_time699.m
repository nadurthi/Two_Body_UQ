%2 Body problem entropy reconstruction
% function MaxEnt2Bodypdfs_time699
% for the 2BP at time step 699 only the pdfs are reconstructed for all 

load('body2simsOTHERS_final')
w=wcut4;
SEX=Xcut4;
 ntt=size(SEX,1);
 nt=length(1:2:ntt);
 
 MaxEntData_cut4=cell(1,12);
 
%   MaxEntData_mc{1}=zeros(nt,3);
%  MaxEntData_mc{2}=zeros(nt,6);
%  MaxEntData_mc{3}=zeros(nt,6);
%  MaxEntData_mc{4}=zeros(nt,6);
%  MaxEntData_mc{5}=zeros(nt,6);
%  MaxEntData_mc{6}=zeros(nt,6);
%  MaxEntData_mc{7}=zeros(nt,10);
%  MaxEntData_mc{8}=zeros(nt,10);
%  MaxEntData_mc{9}=zeros(nt,10);
%  MaxEntData_mc{10}=zeros(nt,15);
%  MaxEntData_mc{11}=zeros(nt,15);
%  MaxEntData_mc{12}=zeros(nt,15); 
 % 1- scaling d xy yz xz
 % 2-min bnd xy yz xz
 % 3- max bnd xy yz xz
 % 4-lam_2xy
 % 5-lam_2yz
 % 6-lam_2xz
 % 7-lam_3xy
 % 8-lam_3yz
 % 9-lam_3xz
 % 10-lam_4xy
 % 11-lam_4yz
 % 12-lam_4xz
 dscale=zeros(nt,6);
 minbnd_xy_yz_xz=zeros(nt,6);
 maxbnd_xy_yz_xz=zeros(nt,6);
 LAM2xy=zeros(nt,6);
 LAM2yz=zeros(nt,6);
 LAM2xz=zeros(nt,6);
 LAM3xy=zeros(nt,10);
 LAM3yz=zeros(nt,10);
 LAM3xz=zeros(nt,10);
 LAM4xy=zeros(nt,15);
 LAM4yz=zeros(nt,15);
 LAM4xz=zeros(nt,15); 
T=1:2:ntt;
 parfor tt=1:length(T)
pp=T(tt);
     % Change here   ********************    %%
DXY=SEX(pp,[1,2],:);
DYZ=SEX(pp,[2,3],:);
DXZ=SEX(pp,[1,3],:);

 %    ********************    %%
 
 
ns=size(DXY,2);
N=size(DXY,3);

XY=zeros(N,ns);
YZ=zeros(N,ns);
XZ=zeros(N,ns);

 for i=1:1:N
 XY(i,1)=DXY(1,1,i);
 XY(i,2)=DXY(1,2,i);
 YZ(i,1)=DYZ(1,1,i);
 YZ(i,2)=DYZ(1,2,i);
 XZ(i,1)=DXZ(1,1,i);
 XZ(i,2)=DXZ(1,2,i);
 end 

bndl_xy=[min(XY(:,1)),min(XY(:,2))];
bndu_xy=[max(XY(:,1)),max(XY(:,2))];

bndl_yz=[min(YZ(:,1)),min(YZ(:,2))];
bndu_yz=[max(YZ(:,1)),max(YZ(:,2))];

bndl_xz=[min(XZ(:,1)),min(XZ(:,2))];
bndu_xz=[max(XZ(:,1)),max(XZ(:,2))];

 [XYt,dxy]=transform_domain(XY,bndl_xy,bndu_xy,-1*ones(1,ns),1*ones(1,ns));
 [YZt,dyz]=transform_domain(YZ,bndl_yz,bndu_yz,-1*ones(1,ns),1*ones(1,ns));
 [XZt,dxz]=transform_domain(XZ,bndl_xz,bndu_xz,-1*ones(1,ns),1*ones(1,ns));
 

 [y1,M1xy]=Cal_moments_samples(XYt,w,1,'raw');
 [y2,M2xy]=Cal_moments_samples(XYt,w,2,'raw');
 [y3,M3xy]=Cal_moments_samples(XYt,w,3,'raw');
 [y4,M4xy]=Cal_moments_samples(XYt,w,4,'raw');

 [y1,M1yz]=Cal_moments_samples(YZt,w,1,'raw');
 [y2,M2yz]=Cal_moments_samples(YZt,w,2,'raw');
 [y3,M3yz]=Cal_moments_samples(YZt,w,3,'raw');
 [y4,M4yz]=Cal_moments_samples(YZt,w,4,'raw');

 [y1,M1xz]=Cal_moments_samples(XZt,w,1,'raw');
 [y2,M2xz]=Cal_moments_samples(XZt,w,2,'raw');
 [y3,M3xz]=Cal_moments_samples(XZt,w,3,'raw');
 [y4,M4xz]=Cal_moments_samples(XZt,w,4,'raw');
 

%%%%%%%%%%%%%%% 2nd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  s=1;
  
y=[zeros(1,ns);y1;y2];

M=[1;M1xy;M2xy];
[y,lam2xy,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));

M=[1;M1yz;M2yz];
[y,lam2yz,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));

M=[1;M1xz;M2xz];
[y,lam2xz,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));
disp('2 moms done')
%%%%%%%%%%%%%%% 3rd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3];
 
M=[1;M1xy;M2xy;M3xy];
[y,lam3xy,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));

M=[1;M1yz;M2yz;M3yz];
[y,lam3yz,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));

M=[1;M1xz;M2xz;M3xz];
[y,lam3xz,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));
disp('3 moms done')
%%%%%%%%%%%%%%% 3rd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4];

M=[1;M1xy;M2xy;M3xy;M4xy];
[y,lam4xy,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));

M=[1;M1yz;M2yz;M3yz;M4yz];
[y,lam4yz,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));

M=[1;M1xz;M2xz;M3xz;M4xz];
[y,lam4xz,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns));
disp('4 moms done')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % 1- scaling d xy yz xz
   
 % 2-min bnd xy yz xz
 % 3- max bnd xy yz xz
 
 % 4-lam_2xy
 % 5-lam_2yz
 % 6-lam_2xz
 
 % 7-lam_3xy
 % 8-lam_3yz
 % 9-lam_3xz
 
 % 10-lam_4xy
 % 11-lam_4yz
 % 12-lam_4xz
 dscale(tt,:)=[dxy,dyz,dxz];
 
 minbnd_xy_yz_xz(tt,:)=[bndl_xy,bndl_yz,bndl_xz];
 maxbnd_xy_yz_xz(tt,:)=[bndu_xy,bndu_yz,bndu_xz];
 
 LAM2xy(tt,:)=lam2xy;
 LAM2yz(tt,:)=lam2yz;
 LAM2xz(tt,:)=lam2xz;
 
 LAM3xy(tt,:)=lam3xy;
 LAM3yz(tt,:)=lam3yz;
 LAM3xz(tt,:)=lam3xz;
 
 LAM4xy(tt,:)=lam4xy;
 LAM4yz(tt,:)=lam4yz;
 LAM4xz(tt,:)=lam4xz; 
 
  SAVE_SHIT(pp,[dxy,dyz,dxz],[bndl_xy,bndl_yz,bndl_xz],[bndu_xy,bndu_yz,bndu_xz],lam2xy,lam2yz,lam2xz,lam3xy,lam3yz,lam3xz,lam4xy,lam4yz,lam4xz) ;
 end
  MaxEntData_cut4{1}=dscale;
 
 MaxEntData_cut4{2}= minbnd_xy_yz_xz;
 MaxEntData_cut4{3}= maxbnd_xy_yz_xz;
 
 MaxEntData_cut4{4}=LAM2xy;
 MaxEntData_cut4{5}=LAM2yz;
 MaxEntData_cut4{6}=LAM2xz;
 
 MaxEntData_cut4{7}=LAM3xy;
 MaxEntData_cut4{8}=LAM3yz;
 MaxEntData_cut4{9}=LAM3xz;
 
 MaxEntData_cut4{10}=LAM4xy;
 MaxEntData_cut4{11}=LAM4yz;
 MaxEntData_cut4{12}=LAM4xz; 
 save('CUT4_MaxEnt_fullData_Dec2012','MaxEntData_cut4','-v7.3')
 matlabpool close
end
%% 
function SAVE_SHIT(tt,dscale,minbnd_xy_yz_xz,maxbnd_xy_yz_xz,LAM2xy,LAM2yz,LAM2xz,LAM3xy,LAM3yz,LAM3xz,LAM4xy,LAM4yz,LAM4xz) 
MaxEntData_cut4=cell(1,12); 
MaxEntData_cut4{1}=dscale;
 
 MaxEntData_cut4{2}= minbnd_xy_yz_xz;
 MaxEntData_cut4{3}= maxbnd_xy_yz_xz;
 
 MaxEntData_cut4{4}=LAM2xy;
 MaxEntData_cut4{5}=LAM2yz;
 MaxEntData_cut4{6}=LAM2xz;
 
 MaxEntData_cut4{7}=LAM3xy;
 MaxEntData_cut4{8}=LAM3yz;
 MaxEntData_cut4{9}=LAM3xz;
 
 MaxEntData_cut4{10}=LAM4xy;
 MaxEntData_cut4{11}=LAM4yz;
 MaxEntData_cut4{12}=LAM4xz; 
save(strcat('CUT4_MaxEnt_partialData_Dec2012_time_',num2str(tt)),'MaxEntData_cut4','-v7.3')
end
  %% x-y plot
%   c=12;
%   
%    [xx,zz]=meshgrid(linspace(-1*s,1*s,200),linspace(-1*s,1*s,200));
%    pent=zeros(size(xx));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         for k=1:1:length(W)
%  pent(i,j)= pent(i,j)+pdf_MaxEnt([xx(i,j),zz(i,j),Xq(k)],lam,y);
%         end
%     end
%  end
% %  c=pdf_MaxEnt(M1,lam,y);
%  figure
%  plot(XYpts(:,1),XYpts(:,2),'ro')
%  xlabel('x')
%  ylabel('z')
%  hold on
%   contour(xx,zz,pent,linspace(0.01,12,20))
%    plot_prop_paper
%  axis([-s,s,-s,s])
%   
%   
%   figure
%    mesh(xx,zz,pent)
%  hold on
%  xlabel('x')
%  ylabel('z')
%   plot_prop_paper
% axis([-s,s,-s,s])
%   
%     %% y-z plot
%    [xx,zz]=meshgrid(linspace(-1*s,1*s,200),linspace(-1*s,1*s,200));
%    pent=zeros(size(xx));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         for k=1:1:length(W)
%  pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([Xq(k),xx(i,j),zz(i,j)],lam,y);
%         end
%     end
%  end
% %  c=pdf_MaxEnt(M1,lam,y);
%  figure
%  plot(XYpts(:,2),XYpts(:,3),'ro')
%  xlabel('y')
%  ylabel('z')
%  hold on
%   contour(xx,zz,pent,linspace(0.01,c/20,20))
%      plot_prop_paper
%  axis([-s,s,-s,s])
%   
%   
%   figure
%    mesh(xx,zz,pent)
%  hold on
%  xlabel('y')
%  ylabel('z')
%    axis([-s,s,-s,s])
%    plot_prop_paper   
%     %% x-z plot
%    [xx,zz]=meshgrid(linspace(-1*s,1*s,200),linspace(-1*s,1*s,200));
%    pent=zeros(size(xx));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         for k=1:1:length(W)
%  pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([xx(i,j),Xq(k),zz(i,j)],lam,y);
%         end
%     end
%  end
% %  c=pdf_MaxEnt(M1,lam,y);
%  figure
%  plot(XYpts(:,1),XYpts(:,3),'ro')
%  xlabel('x')
%  ylabel('z')
%  hold on
%   contour(xx,zz,pent,linspace(0.00001,c/5,20))
%      plot_prop_paper
%  axis([-s,s,-s,s])
%   
%   
%   figure
%   mesh(xx,zz,pent)
%  hold on
%  xlabel('x')
%  ylabel('z')
%    axis([-s,s,-s,s,0,2*c])
%       plot_prop_paper