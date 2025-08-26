%% Author: Enping Lin, Ze Fang
%% Contact: enpinglin@xmu.edu.com
%% Date: 2025-08-26
%% This script is to reproduce 1H-1H symmetrical spectra reconstruction using IST and SCREEN.
%% The spectra data are obtained from the 128×128 COSY projection of the 3D TOCSY-COSY 
%% spectrum of quinine.
%%
clear all,clc, close all
addpath('l1_ls_matlab');
addpath('data')
addpath('function')
load('quinine.mat')

%% Read data
Fid = ifft2(Spec2D_Ideal);
[N1,N2] = size(Fid);
t1 = 0:N1-1;
t2 = 0:N2-1;
ppm1 = linspace(-1.6,8.37,N1);
ppm2 = linspace(-1.6,8.37,N2);
%% ---------- generate mask -----------------------------
ist_screen = 1;   % 0 for SCREEN, 1 for IST
ty = 2;           % 1 for random sampling schedule, 2 for 2D poisson sampling schedule,
                  % 3 for SCPG sampling schedule
if ty == 1
    load mask_ran.mat
elseif ty == 2
    load mask_2Dpoisson.mat
elseif ty == 3
    load mask_sympoisson.mat
end

maskc = mask;
FidNus_fill = Fid;
%--------- Symmetrically Copy ----------------------------- 
if ty == 3
    maskc = logical(mask'+mask); 
    FidNus = Fid.*mask;
    FidNus_fill = (FidNus+FidNus.').*(ones(N1)-0.5*diag(ones(1,N1)));
end
%%
FidNus_fill(maskc==0) = [];
FidNus_fill = FidNus_fill(:);
sr = length(find(mask==1))/numel(mask);

%%
tic
options.max_iter = 10000;
options.res = 1e-5;
options.change_thr = 1;
%%
if ist_screen == 1
    X_Output = ist_d2D_LEP(FidNus_fill, maskc, options); 
else
    ye = [real(FidNus_fill);imag(FidNus_fill)];
    NorFac = max(abs(ye));
    ye = ye./NorFac;
    AM = Fu_downsample(maskc,N1,N2);
    AMh = AM';
    [SpecRec,status,history] = l1_ls(AM,AMh,length(ye),N1*N2*2,ye,0.001,1e-5); % 1e-5
    flat_spectrum = SpecRec(1:length(SpecRec)/2)+1i*SpecRec(length(SpecRec)/2+1:end);
    flat_spectrum = flat_spectrum * NorFac;
    X_Output = reshape(flat_spectrum,N1,N2);
end 

toc
X_Output = real(X_Output);
X_Output = X_Output./max(X_Output(:));
%% show 2D result
TyC = [0.47,0.67,0.19;0.49,0.18,0.56;0.00,0.00,1.00;  0.85,0.33,0.10;0.64,0.08,0.18];
FSC = [0,0.4,0.7]; 

Spec2D_Ideal = Spec2D_Ideal./max(real(Spec2D_Ideal(:)));
figure,contour(ppm1,ppm2,(Spec2D_Ideal),linspace(0.03,1,50),'LineColor',FSC ), %title('Full Sampling')
set (gca, 'XColor',FSC ,'YColor',FSC ,'box','on','linewidth',1.5,'fontsize',20,'fontname','Calibri','fontweight','bold','xdir','reverse','ydir','reverse');
figure,contour(ppm1,ppm2,X_Output,linspace(0.03,1,50),'LineColor',TyC(ty,:));%title('ist NUS Reconstruction ')
set (gca,'XColor',TyC(ty,:),'YColor',TyC(ty,:),'box','on','linewidth',1.5,'fontsize',20,'fontname','Calibri','fontweight','bold','xdir','reverse','ydir','reverse');

%% Diagonal Peak analysis
RecS = X_Output;
IdeS = Spec2D_Ideal;
[pks,locs,w,p] = findpeaks(diag(RecS),'minpeakheight',0.03);
bnd1 = max(floor(locs-w/2),1);
bnd2 = min(ceil(locs+w/2),N1);
IdeS_diagPInt = zeros(1,length(locs));
RecS_diagPInt = zeros(1,length(locs)); 
for it = 1:length(locs)
    IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
    RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
end
CorDiagPeaks= corrcoef(IdeS_diagPInt, RecS_diagPInt);
Dia_RLNE = norm(IdeS_diagPInt-RecS_diagPInt)/norm(IdeS_diagPInt);
sl = max([IdeS_diagPInt,RecS_diagPInt])*1.1;
Dia_lx = polyfit(IdeS_diagPInt,RecS_diagPInt,1);
figure,plot(IdeS_diagPInt,RecS_diagPInt,'o','MarkerFaceColor','b','MarkerEdgeColor','k'), xlim([0,sl]), ylim([0,sl])
hold on, line([0,sl],[0*Dia_lx(1)+Dia_lx(2),sl*Dia_lx(1)+Dia_lx(2)],'LineWidth',2,'color','b')
set(gca,'fontsize',18,'fontweight','bold','fontname','times new romance','LineWidth',2) 

%% ================ find cross peak Index =========================
[ B ] = BDiagOnes( N1,10 );
DiagFlag = logical(  B  );
IdeS_RevDiag = IdeS.*(~DiagFlag);
RecS_RevDiag = RecS.*(~DiagFlag);

bwm_Ideal = ((IdeS_RevDiag+IdeS_RevDiag.')/2>0.03);
bwm_Ideal = bwareaopen(bwm_Ideal,1) ;
[L_Ideal,num_Ideal] = bwlabel(bwm_Ideal);
IdePeak = zeros(1,length(num_Ideal));
RecPeak = zeros(1,length(num_Ideal));
for it = 1: num_Ideal
    Temp_Ideal = find(L_Ideal==it);
    IdePeak(it) = sum(IdeS(Temp_Ideal ),'all');
    RecPeak(it) = sum(RecS(Temp_Ideal ),'all');
end    

CorCroPeaks = corrcoef( IdePeak, RecPeak);
Cro_lx = polyfit(IdePeak,RecPeak,1);
Cro_RLNE = norm(IdePeak-RecPeak)/norm(IdePeak);
lmax = max([IdePeak,RecPeak]);
fl = lmax*1.1;

figure,subplot(2,1,1)
stem(sort(IdePeak),'marker','_','linewidth',2,'Color',FSC ),
set(gca,'fontsize',18,'fontweight','bold','fontname','times new romance','LineWidth',2)
xlim([0,length(RecPeak)+1]),ylim([0,fl])
subplot(2,1,2)
stem(sort(RecPeak),'marker','_','linewidth',2,'Color',TyC(ty,:)),
set(gca,'fontsize',18,'fontweight','bold','fontname','times new romance','LineWidth',2)
xlim([0,length(RecPeak)+1]),ylim([0,fl])

figure,plot( IdePeak,RecPeak,'o','MarkerFaceColor','r','MarkerEdgeColor','k'), xlim([0,fl]), ylim([0,fl])
hold on, line([0,fl],[0*Cro_lx(1)+Cro_lx(2),fl*Cro_lx(1)+Cro_lx(2)],'LineWidth',2,'color','r')
set(gca,'fontsize',18,'fontweight','bold','fontname','times new romance',...
    'LineWidth',2) 
title('CrossPeaks Intensity')
