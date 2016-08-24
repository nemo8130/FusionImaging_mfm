% FusionImaging_mfm is a MATLAB code for processing and analyses of Images 
% generated from Magnetic Force Microscopy (MFM) 
% The code was compiled in MATLAB Version: 9.0.0.341360 (R2016a) 
% Users are recommended to use this or higher version with the 
% 'Image Processing Toolbox' installed.
%
% The software converges a collection of raw MFM images for a given sample 
% collected at different Lift-heights into Fusion Images. 
% 
% The function needs two and exactly two input arguments
% Usage : FusionImaging_mfm (sampfold, datmm)
% where 
% sampfold: Is the sample folder with full path containing raw MFM .tif images\n')
% and 
% datmm: Is the data file containing Minimum and Maximum voltages (in units of mV) 
% of each image in the 'sampfold' stored in the tabular format specified below. 
%
% Directory            Filename                     Min_Voltage(mV)         Max_Voltage(mV)
% APTES_1_P             ap_010.tif                   -1200.00                -1100.00
% APTES_1_P             ap_020.tif                   -1100.00                  979.80
% APTES_1_P             ap_030.tif                    -970.20                  875.00
% APTES_1_P             ap_040.tif                    -821.60                  763.50
% APTES_1_P             ap_050.tif                    -821.60                  763.50
% APTES_1_P             ap_060.tif                    -730.40                  681.00
% APTES_1_P             ap_070.tif                    -824.30                  750.50
% BSA_1_P               bs_060.tif                      -4.90                    4.90
% BSA_1_P               bs_080.tif                      -4.40                    4.50
% BSA_1_P               bs_090.tif                      -4.70                    4.60
% CITRATE_1_P           c1_010.tif                   -3200.00                 2500.00
% CITRATE_1_P           c1_020.tif                   -4600.00                 4000.00
% CITRATE_1_P           c1_030.tif                   -1900.00                 1500.00
% CITRATE_1_P           c1_040.tif                   -3400.00                 2600.00
% CITRATE_1_P           c1_050.tif                   -3300.00                 2400.00
% CITRATE_1_P           c1_060.tif                   -3700.00                 2600.00
% CITRATE_2_P           c2_010.tif                    -174.70                  135.30
% CITRATE_2_P           c2_020.tif                    -163.90                  129.70
% CITRATE_2_P           c2_030.tif                    -169.70                  126.00
%
% An example datmm file (Min_Max_Table_example_phase_data.format) is provided with this distribution.
% Note: The datmm file should contain the commented header starting with a '%'
% Images must have '.tif' extension
%
% Specification of Full path needs to be consistent with the OS (Unix / Windows: see example 'sampfold')
% 
% Full path may or may not end with a '/' (Linux) or '\' (Windows) : The program can handle both
% the datmm file can be directly called from the current directory by the file name alone
%
% logz=FusionImaging_mfm (sampfold, datmm) will return the (NxN) matrix containing 
% coefficients of the inverse-square term of the fitted polynomial; 
% (magnitudes scaled to natural logarithm, log_e)
% These coefficients are representative of the magnetic fields 
% corresponding to each small square grid in the 'two dimensional' image-space. 
%
% Example Input:  
%
% sampfold='/home/sankar/Dropbox/sankar-puja/Data/Data_Widout_VW/phase_data/' (Linux / Unix) 
% sampfold='C:\user\Sankar\Desktop\mfmdata\Data_Widout_VW\phase_data\'; (Windows / PC)
%
% datmm='/home/sankar/Dropbox/sankar-puja/Data/SOFTWARE/Min_Max_Table_example_phase_data.format'; (File with full path: Linux)
% datmm='Min_Max_Table_example_phase_data.format';  (File kept in the Current Directory)
%
% datmm='C:\user\Sankar\Data\SOFTWARE\Min_Max_Table_example_phase_data.format'; (File with full path: Windows)
%
%

function logz=FusionImaging_mfm(sampfold,datmm)
blk=2;      % Image_Dimension: 400*400*1*1 = 160000 sq_unit = 2*2*200*200
A={};B={};H={};G={};N={};NW={};
subplot(2,2,1)
k=1;       % Only run for the chosen input

% Identify the OS (Linux / Windows)

if (isunix == 1 & ispc == 0)           % Linux 
    g1=strfind(sampfold,'/');
elseif (isunix == 0 & ispc == 1)       % Windows
    g1=strfind(sampfold,'\');
end

lg1=length(g1)

% Identify whether the string ends with an oblic ('/': Linux; '\': Windows)
% And act accordingly

if (length(sampfold) == g1(lg1))    
    dinp=sampfold(g1(lg1-1)+1:g1(lg1)-1);
    %fprintf('I am in if\n')
elseif (length(sampfold) > g1(lg1))
    dinp=sampfold(g1(lg1)+1:length(sampfold));    
    if (isunix == 1 & ispc == 0)                    % Linux 
        sampfold=strcat(sampfold,'/');
    elseif (isunix == 0 & ispc == 1)                % Windows
        sampfold=strcat(sampfold,'\');
    end        
    %fprintf('I am in elseif\n')
end

dinp
sampfold

[a,b]=imxtract(sampfold,'tif');
A{k}=a;B{k}=b;
h=liftheights(a);
H{k}=h';

%fprintf('%10d  %15s\n',i,dinp);
[newimA,newimB,newim,txt,h1,X,N,phi]=mfmanal_readfile_sankar(b,a,blk,dinp,datmm);
NW{k}=newim;

ppp=[];

for ii=1:length(NW{k}{1})
    [ppN,x,ybig]=bigP(NW,H,k,ii);
    ppp=[ppp;ppN(1)];
end

rsz=sqrt(length(ppp));
csz=rsz;

z=reshape(ppp,[rsz csz]);
zmax=max(ppp);
zmin=min(ppp);
dz1=zmin+(zmax-zmin)/2;

%fprintf('Here I am: So the control statement have not worked\n');

z1I=zeros(size(z));
z2I=z1I;
z3I=z1I;

z1I(find(z>1.5*dz1))=z1I(find(z>1.5*dz1))+1;        % Red matrix
z2I(find((z>=0.75*dz1) | (z<=1.5*dz1)))= z1I(find((z>=0.75*dz1) | (z<=1.5*dz1))) + 1;    % Green Matrix
z3I(find(z<0.75*dz1))=z3I(find(z<0.75*dz1))+1;      % Blue matrix

RGBI = [];

RGBI(:,:,1) = z1I;
RGBI(:,:,2) = z2I;
RGBI(:,:,3) = z3I;

subplot(2,2,2)

RI = imref2d(size(RGBI));
imshow(RGBI,RI)
title (['one-zero 2D:  ',num2str(strrep(dinp,'_','-'))])
set(gca,'fontsize',15)
colorbar

%========================= Added part by Sankar ======================

z1=zeros(size(z));
z2=zeros(size(z));
z3=zeros(size(z));

[r1 c1] = find(z>1.5*dz1);
[r2 c2] = find((z>=0.75*dz1) | (z<=1.5*dz1));            % Check and Invert
[r3 c3] = find(z<0.75*dz1);

if (length(r1)>0)
    for i = 1:length(r1)
        z1(r1(i),c1(i))=abs(z(r1(i),c1(i)));           % abs ?
        %fprintf('%10.2f%2s%10.2f\n',z(r1(i),c1(i)),'->',z1(r1(i),c1(i)));
    end
end

if (length(r2)>0)
    for i = 1:length(r2)
        z2(r2(i),c2(i))=abs(z(r2(i),c2(i)));           % abs ?
        %fprintf('%10.2f%2s%10.2f\n',z(r1(i),c1(i)),'->',z1(r1(i),c1(i)));
    end
end

if (length(r3)>0)
    for i = 1:length(r3)
        z3(r3(i),c3(i))=abs(z(r3(i),c3(i)));           % abs ?
        %fprintf('%10.2f%2s%10.2f\n',z(r3(i),c3(i)),'->',z3(r3(i),c3(i)));
    end
end

z1max = max(max(abs(z1)));
z2max = max(max(abs(z2)));
z3max = max(max(abs(z3)));

if (z1max > 0)
    z1 = z1./z1max;
end

if (z2max > 0)
    z2 = z2./z2max;
end

if (z3max > 0)
    z3 = z3./z3max;
end

z1=uint8(round(z1*255));
z2=uint8(round(z2*255));
z3=uint8(round(z3*255));

RGB = [];

RGB(:,:,1) = z1;
RGB(:,:,2) = z2;
RGB(:,:,3) = z3;

subplot(2,2,3)
RI = imref2d(size(RGB));
imshow(RGB,RI)
title (['RGB 2D:  ',num2str(strrep(dinp,'_','-'))])
set(gca,'fontsize',15)
colorbar

subplot(2,2,4)
logz = log(abs(z));
imagesc(logz)
axis equal
axis ([0 200 0 200])
ax=gca;
ax.XTickLabel = {'0','20','40','60','80','100','120','140','160','180','200'};
ax.YTickLabel = ax.XTickLabel;
title (['log(z) 2D:  ',num2str(strrep(dinp,'_','-'))])
set(gca,'fontsize',15)
colorbar

Y = 1;
end

%=================== Subroutines begins here =============================

function SubB=im2blk(A,blks)
if ischar(A)==1
    A=imread(A);
end
% i assume 8-bit gray scale image
[m,n,k]=size(A); % and m=n with 1 channel k=1
ImageSize=m*n;
BlockD=blks;
BlockSize=BlockD*BlockD;
NoOfBlock=ImageSize/BlockSize;
SubB=zeros(BlockD,BlockD,NoOfBlock); %arrays of blocks.
B=double(A); %important to convert uint8 to double when dialing with image.
% thats what ru asking for.
k=1;
for i=1:BlockD:m
    for j=1:BlockD:n
        SubB(:,:,k)=B(i:i+BlockD-1,j:j+BlockD-1); k=k+1;
    end
end
%compare between first submatrix A with first block.. its the same elements.
%B(1:2,1:2)
%SubB(:,:,1)
end


function [pp,x,ybig]=bigP(NW,H,k,l);
w=NW{k};
x=H{k};
ybig(1)=0;
for i=1:length(x)
    y=w{i};
    s=y(:,:,l);
    ybig(i)=mean(mean(s));
    
end
x=x(:);
ybig=ybig(:);

xx=linspace(min(x),max(x),100);
if (length(x)<4)
    [P,s]=polyfit(x,ybig,1);
    yy=P(1)*xx+P(2);
else
    [P,s]=polyfit(x,ybig,3);
    yy=P(1)*xx.^3+P(2)*xx.^2+P(3)*xx+P(4);
end
xxx=1./xx;
[pp,ss]=polyfit(xxx,yy,2);       % Try out '4'
end


function [pp,x,ybig]=bigPplot(NW,H,k,l,dinp);
w=NW{k};
x=H{k};
ybig(1)=0;
for i=1:length(x)
    y=w{i};
    s=y(:,:,l);
    ybig(i)=mean(mean(s));
    
end
x=x(:);
ybig=ybig(:);

xx=linspace(min(x),max(x),100);
if (length(x)<4)
    [P,s]=polyfit(x,ybig,1);
    yy=P(1)*xx+P(2);
else
    [P,s]=polyfit(x,ybig,3);
    yy=P(1)*xx.^3+P(2)*xx.^2+P(3)*xx+P(4);
end
xxx=1./xx;
[pp,ss]=polyfit(xxx,yy,2);

%length(ybig)

if (length(ybig) > 3)                              % At least 4 points required to perform the regstats
    stats = regstats(xxx,yy,'quadratic','rsquare')
    Rsq = stats.rsquare;
else
    Rsq = 0.005;
end

plot(x,ybig,'ro',xx,yy,'b','MarkerSize',5,'LineWidth',3.5);
xlabel 'Lift-Heights'
ylabel 'Force-Equivalent'
P1form=sprintf('%20.10f\n',P(1));
title (['Polyfit:  ',num2str(strrep(dinp,'_','-')),';  P1=',num2str(P(1)),'; Rsq=',num2str(Rsq)])    % If you want a float format, replace P(1) by P1form
set(gca,'fontsize',20)
end


function [newimA,newimB,newim,txt,h,X,N,phi,flag]=mfmanal_readfile_sankar(image_cell,image_cellname,blk,dinp,datmm);
% Image cell is a cell contains image stacks of different heights h
flag=0;
s=image_cell;
h=liftheights(image_cellname);
h=h(:);
newimA={};newimB={};X={};newim={};N={};
l=length(s);

fprintf('%18s  %6d\n','Number of images:',l);

for i=1:l
    z=s{i};
    finp=image_cellname{i};
    fprintf('%20s  %20s  %20s  %20s\n','image_directory:',dinp,'image_filename:',finp);
    if length(size(z))>2
        z=rgb2gray(z);
    end
    [r,c]=size(z);
    newimA{i}=z(30:429,30:429);
    G=z(1:r,r-10:c);
    newimB{i}=G;
    [n,x]=imhist(newimA{i});
    nf = n./sum(n);
    plot(x,nf)
    AUC=trapz(x,nf);
    %fprintf('Trapizoid Area Under the Imhist Relative Intensity vs. Pixel Curve: %10.5f\n',AUC);
    xlabel 'Scaled Pixel Values'
    ylabel 'Relative Intensity'
    title (['Polydispersity: ',num2str(strrep(dinp,'_','-'))])
    set(gca,'fontsize',15)
    if i==1;hold on;end
    N{i}=n;
    %     X{i}=x;
    ik=find(n==max(n));
    if length(ik)>1
        x=mean(x(ik));
    else
        x=x(ik);
    end
    X{i}=x;
    
    [minV,maxV]=returnMinMax(dinp,finp,datmm);
    
    txt{i}=[minV maxV];
    X{i}=X{i}*((maxV-minV)/255)+minV;
    G=double(newimA{i});
    G=G*((maxV-minV)/255)+minV;
    newim{i}=im2blk(G,blk);
end
phi=[];
for i=1:length(X)
    phi=[phi;X{i}];
end
end


function [fname,cw]=imxtract(path,extension)
fname={};cw={};
z=dir([path '*.' extension ]);
frac=[];
k=length(z);
count=1;
for i=1:k
    fname{i}=z(i,1).name;
    cw{i}=imread([path fname{i}]);
end
end


function h=liftheights(b)
% sampfold is the directory
% b contains all the file names
l=length(b);
h=[];
for i=1:l
    s=b{i};
    g=strfind(s,'_');
    l=strfind(s,'.');
    k=s(g+1:l-1);
    k=str2num(k);
    h=[h;k];
end
end

function [minVmv,maxVmv]=returnMinMax(dinp,finp,datmm);

%datmm='Min_Max_Table_example_phase_data.format';

minmax=importdata(datmm);
minmaxV=minmax.data;
rows=minmax.textdata;
header=rows(1,:);
dirname=rows(2:length(rows),1);
fname=rows(2:length(rows),2);

minV=minmaxV(:,1);
maxV=minmaxV(:,2);

%whos dirname fname minV maxV

% dinp=input('Enter dirname:\n','s');
% finp=input('Enter filename:\n','s');

flag=0;
for i = 1:length(dirname)
    dirstr=dirname{i};
    fnamestr=fname{i};
    %whos str1 dinp
    match1=strcmp(dirstr,dinp);
    match2=strcmp(fnamestr,finp);
    %disp(match1);disp(match2)
    flag=0;
    if (match1 == 1 & match2 == 1)
        minVmv=minV(i);
        maxVmv=maxV(i);
        %fprintf('%20s  %5d  %15s  %15s  %10.3f  %10.3f\n','Match found at i=',i,dirstr,fnamestr,minVmv,maxVmv)
        flag=1;
        break;
    end
end

    if (flag==0)
    fprintf('Data not found for %20s %20s\n',dirstr,fnamestr)
    return
    end
end



