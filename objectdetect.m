function [objN objco Igrm AgN0 els2 perim sED2 dmaxO sSOL2 subimg Graygrain_whole_selected] = OBJECTDETECT(file,glevel,rK,Aco);
%% OBJECTDETECT: A program to detect objects from a Force-Microscopy Image
% And filter objects based on the set of criteria enlisted below
%
% 1. Thresholding: A certain 'level' is used to set the grayscale threshold
% which is then used to delineate the boundary between the objects (white) and the background (black). 
% Unless specifically guided by the user, the program determines this level from the 'graythresh' function 
% applied on the particular input image.
%
% 2. Shape : A cutoff given on the roundedness of each object to discard unusually asymmetric ones.
% Pseudo objects or image artefacts (scratches for example) are elegantly
% dealt with this filter.
%
% 3. Area : Small (dotted) objects are discarded based on a filter on the raw area
% (in pixel^2) 
%
%-----------------------------------------------------------------------------------------------------
%% USAGE: 
% [objN objco Igrm AgN0 els2 perim sED2 dmaxO sSOL2 subimg Graygrain_whole_selected]=OBJECTDETECT(file,glevel,rK,Aco);
% 
% where 
% The Input Arguments mean the following:
%
% file: The Image filename with full path or just the filename for a file residing in the current directory
%
% glevel: The threshold level (suggested value: 0.40); A choice of '0' for this parameter will restore the default
% i.e., the value for parameter will be returned from graythresh
%
% rK: The tuning parameter for roundedness: The program fits each object to an ellipse 
% and calculates the ratio (r) of major-to-minor-axis-lengths. 
% From the distribution of the r values it then calculates the mean (mu(r)) and the standard deviation
% (sigma(r)) and sets a cut-off on the higher end of the distribution to filter out objects based on shape. 
% i.e., only those objects are retained, for which, (1 <= r <= (mu(r) + rK*sigma(r))). 
% (Default for rK: 1.5) [suggested range: 1.0 to 2.0]
% To note, the lower bound of r is 1 by definition.
%
% Aco: is the area cut-off in pixel^2 to filter out small dotted objects (Default: 50)
%---------------------------------------------------------------------------------------
% Return variables:
% objN:  the array containing the object numbers
% objco: coordinates (x,y) of the object-centroid.
% Igrm:  median of the unscaled pixel intensities corresponding to the object
% AgN0:  Area (in pixel^2) 
% els2:  Roundedness 
% perim: Perimeter 
% sED2:  Equivalent Diameter 
% dmax0: Maximum Diameter
% sSOL2: Solidity 
% subimg: cell array containing the subarrays corresponding to each
% filtered object 
% Graygrain_whole_selected: All selected objects (for visualisation)
%--------------------------------------------------------------------------------------
% To view the i'th filtered object on the background of the original 'full' image, 
% run: IMSHOW(subimg{i})
% and render it using the coordinates of the object centroids:
%
% hold on;plot(objco(i,1),objco(i,2),'r.','MarkerSize',5)
% hold on;plot(objco(i,1),objco(i,2),'yo','MarkerSize',50)
%
% or view them one by one with pseudo color
% figure;for i = 1:length(objN);IMAGESC(subimg{i});pause(1);end
% or in grayscale
% figure;for i = 1:length(objN);IMAGESC(subimg{i});colormap('gray');pause(1);end
% 
% or view all selected objects
% 
% figure; imagesc(Graygrain_whole_selected);colormap('gray'); axis equal; axis([0.5 518 0 500])
%----------------------------------------------------------------------------------------------------
%% EXAMPLE: 
% (with suggested default values for the inputs)
%
% [objN objco Igrm AgN0 els2 perim sED2 dmaxO sSOL2 subimg Graygrain_whole_selected]=OBJECTDETECT('adg-testim-addborder.png',0,1.5,50);
% 
%---------------------------------------------------------------------------------------------------

outfile = strcat(file,'.obj')
fid1=fopen(outfile,'w');
bg = 15;
perrm=0.05;
I = imread(file);
I1 = rgb2gray(I);
Npix = size(I1)     % Assigned as [y x]

%==========================================================================
%========== Identify OS and extract file-tag for title ====================
%==========================================================================

if (isunix == 1 & ispc == 0)           % Linux 
    g1=strfind(file,'/');
elseif (isunix == 0 & ispc == 1)       % Windows
    g1=strfind(file,'\');
end

lg1=length(g1);

% Identify whether the string ends with an oblic ('/': Linux; '\': Windows)
% And act accordingly

if (lg1==0)
    filetag=file;
elseif (length(file) == g1(lg1))
    filetag=file(g1(lg1-1)+1:g1(lg1)-1);
    %fprintf('I am in if\n')
elseif (length(file) > g1(lg1))
    filetag=file(g1(lg1)+1:length(file));
end

ftag=strsplit(filetag,'.');
filetag=char(ftag(1));
filetag=strrep(filetag,'_','-');


%==========================================================================

background = imopen(I1,strel('disk',bg));     % default: 15

%==========================================================================
figure(1)
subplot(2,2,1)
imshow(I)
subplot(2,2,2)
imshow(I1)
subplot(2,2,3)
imshow(background,[])
I2 = I1 - background;
subplot(2,2,4)
imshow(I2)
%==========================================================================

%==========================================================================
figure(2)
subplot(1,2,2)
surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
colorbar
set(gca,'ydir','reverse');
I3 = imadjust(I2);
% I3=I2;
subplot(1,2,1)
imshow(I3);
figure(3)

%==========================================================================
if (glevel == 0)
    level = graythresh(I3);     % Thresholding
else
    level = glevel;             % User defiend
end
%==========================================================================

f1=50;                       % less the value, more the noise (ideal range: 25 to 50)
bw = im2bw(I3,level);
bw = bwareaopen(bw,f1);
subplot(1,2,1)
imshow(bw)
xlabel 'White on Black'
set(gca,'fontsize',20)
cc = bwconncomp(bw,4);
Nobj=cc.NumObjects;
grain = false(size(bw));
grain(cc.PixelIdxList{Nobj}) = true;     % default: 50
labeled = labelmatrix(cc);
RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
subplot(1,2,2)
imshow(RGB_label)
xlabel 'RGB on White'
set(gca,'fontsize',20)
%mtit(['Object Distribution: ',num2str(tag)],'FontSize',20)
%==========================================================================

%%=======================================================================
% Extract Objects (returned as boxes by Default)
%%=======================================================================

graindata = regionprops(cc,'basic');
grain_areas = [graindata([1:Nobj]).Area];
grain_cent=[graindata([1:Nobj]).Centroid];
grain_bb=[graindata([1:Nobj]).BoundingBox];

lgc0=length(grain_areas);
lgc1=length(grain_cent);
lgc2=length(grain_bb);

%==========================================================================
figure(4)
imshow(RGB_label)
hold on
CentObj = [];
BoxObj = [];

for i = 1:2:length([grain_cent])
    xc=grain_cent(i);
    yc=grain_cent(i+1);
    plot(xc,yc,'r*')    
    hold on
    CentObj = [CentObj;xc yc];
end

for i = 1:4:length([grain_bb])
    x1=grain_bb(i);
    y1=grain_bb(i+1);
    h1=grain_bb(i+2);
    w1=grain_bb(i+3);
    rectangle('position',[x1,y1,h1,w1]) 
    hold on
    BoxObj = [BoxObj;x1 y1 h1 w1];
end

len_gr=length([grain_cent]);

%whos CentObj
%whos BoxObj
%==========================================================================
%========================================================================
% Test for a single object
%========================================================================
%

[max_area,idx] = max(grain_areas);               % identify unique objects
grain = false(size(bw));
grain(cc.PixelIdxList{idx}) = true;
cont = bwperim(grain,8);
[yy1,xx1]=find(grain==1);
szI=size(I);
RGBgrain=zeros(szI(1),szI(2),3);
Graygrain=zeros(szI(1),szI(2),1);

    for i = 1:length(xx1)
        RGBgrain(yy1(i),xx1(i),:) = I(yy1(i),xx1(i),:);
        Graygrain(yy1(i),xx1(i),:) = I1(yy1(i),xx1(i),:);
    end

Ag0=bwarea(grain);             % Is the correct one

cont2=imdilate(cont,strel('disk',5));
Ag1=bwarea(cont);
Ag2=bwarea(cont2);
%fprintf('%12.3f %12.3f %12.3f\n',Ag0,Ag1,Ag2)

%==========================================================================
figure(5)
AgN = [];
%whos bw
elstore = [];
%subplot(2,2,1)
%imshow(I)
%subplot(2,2,2)
%imshow(I1)
%subplot(2,2,1)
hold on
szI=size(I);
RGBgrain_whole=zeros(szI(1),szI(2),3);
Graygrain_whole=zeros(szI(1),szI(2),1);
%RGBgrain_store = [];

%======================== Loop over all objects ===========================

eflag=0;

for idx = 1:Nobj
    grain = false(size(bw));
    grain(cc.PixelIdxList{idx}) = true;
    %imshow(grainN);
    contN = bwperim(grain,8);
    %imshow(contN)
    [yy,xx]=find(contN==1);
    [yy1,xx1]=find(grain==1);        
    RGBgrain=zeros(szI(1),szI(2),3);
    Graygrain=zeros(szI(1),szI(2),1);
        for ii = 1:length(xx1)       
            RGBgrain_whole(yy1(ii),xx1(ii),:) = I(yy1(ii),xx1(ii),:);
            RGBgrain(yy1(ii),xx1(ii),:) = I(yy1(ii),xx1(ii),:);       
            Graygrain_whole(yy1(ii),xx1(ii),:) = I1(yy1(ii),xx1(ii),:);
            Graygrain(yy1(ii),xx1(ii),:) = I1(yy1(ii),xx1(ii),:);       
        end
    
    %RGBgrain_store = [RGBgrain_store;RGBgrain];
    % Fit an elepse to the points: Calculate the ratio of the major to
    % minor axis (ecen) and set a filter to filter out spourious
    % objects

    xmin=min(xx);
    ymin=min(yy);
    xmax=max(xx);
    ymax=max(yy);
    
    %fprintf('%5d %5d %5d %5d %5d\n',idx,xmin,xmax,ymin,ymax)
    
    % Calculate a distance matrix for points pertaining to each object 
    % and put a filter on it
    % to be taken into consideration for an ellipse fit
 
    dobj = [];
    
    for i10 = 1:length(xx1)
        for j10 = 1:length(yy1)
        dobj(i10,j10)=sqrt((xx1(i10)-xx1(j10))^2 + (yy1(i10)-yy1(j10))^2);
        end
    end
        
    dmax = max(max(dobj));
%    fprintf('Max length and spread (std) of the the object %d : %f %f %f\n',idx,dmax,std(xx1),std(yy1));
    halflen=(Npix(1)+Npix(2))/(2*2);
    
    if (dmax <= halflen)                   % The maximum distance between two points in an object should not exceed half-length of the image
        if (std(xx1) > 0.5 && std(yy1) > 0.5) % filter out artifactual scattered isolated points forming pseudo objects
            if (length(xx)>=6) 
                plot(xx,yy,'.')
                plot(xx1,yy1,'.')
                ellp=fit_ellipse(xx,yy,1);
                est=ellp.status;
                if (length(est) == 0)
                    %fprintf('ellipse found\n')
                    mjax=ellp.long_axis;
                    mnax=ellp.short_axis;
                    ecen=mjax/mnax;
                    elstore(idx) = ecen;
                    text(ellp.X0_in,ellp.Y0_in,[num2str(idx)])
                    eflag=1;
                else
                    elstore(idx) = 0.00;
                    fprintf('ellipse could not be fit\n')
                end
                
            else
                elstore(idx) = 0.00;
            end
        else
            %fprintf('%s %d %s\n','ellipse-fit will be skipped for object:',idx,'for being severely assymetric (dumble-shaped)')
            %fprintf('%s\n','It is a Pseudo-object')
        end
        
        if (eflag==1)
            AgN=[AgN;bwarea(grain)];
        end
    end
    
end
axis([-25 525 -25 525])
title ('All Objects after Thresholding')

%whos I
%whos RGBgrain
%whos RGBgrain_store

%=============================================

length(elstore);
me=median(elstore(find(elstore>0)));
se=std(elstore(find(elstore>0)));
se=1.0;

%================= Inverse Ecentricity Cutoff =======================
%K=1.25;    % Should be somewhere between 1 to 2
K=rK;       % User defined
%====================================================================

mu2sigp = me+(K*abs(se));
mu2sign = me-(K*abs(se));

[nn00,cc00]=hist(elstore);
snn00=nn00./sum(nn00);

% Ratio of Major-to-Minor Axis can not be less than one (by definition)
% So, a meaningful cut off should only be on the higher end (<=) 

indel=find(elstore<=mu2sigp & elstore>=1);
elsfilt = elstore(indel);
[nn1,cc1]=hist(elsfilt);
snn1=nn1./sum(nn1);

figure(6)
subplot(1,2,1)
[nn0,cc0]=hist(elstore);
snn0=nn0./sum(nn0);
bar(cc0,snn0)
xlabel 'Major-to-Minor-axis-ratio'
ylabel 'Norm Freq'
title 'Distribution of Size of the objects'
subplot(1,2,2)
bar(cc1,snn1)
xlabel 'Major-to-Minor-axis-ratio'
ylabel 'Norm Freq'
title 'Distribution of Size of the objects (after filtering based on size)'

%return

%=========== Area Cutoff on Pixel^2 =====================

[n0,c0]=hist(AgN,20);
sn0=n0./sum(n0);
meA=mean(AgN);
seA=std(AgN);
m2snA=meA-2*seA;

%======================= Area Cutoff (mic^2) =============
%m2snA=0.015;
%=========================================================

figure(7)
bar(c0,sn0,'yellow')
xlabel 'Area (pixel^2)'
ylabel 'Norm Freq'
title 'Area Histogram (before filtering small objects)'

%m2snA=input('Have a look at the Area Histogram \nNow Enter the Area Cutoff (in pixel^2) for filtering small objects:\n');
m2snA=50;         % Optimized on real data
m2snA=Aco;        % User defined
hold on;
bar(m2snA,max(sn0),'k','BarWidth',1.0)
ia=find(c0>m2snA);
fracA=sum(sn0(ia))*100;

fprintf('---------------------------------------------------------------------------\n');
fprintf('%s %6.2f %s %6.2f %s\n','The cutoff given in area (',Aco,' pixel^2) filtered ',fracA,'% of objects discarding small dotted ones');
fprintf('---------------------------------------------------------------------------\n\n');

%=========================================================

%=========================================================================
% Now filter the objects based on (i) shape and (ii) area
%=========================================================================

figure(8)
%subplot(2,2,4)
AgN0 = [];
hold on
fprintf(fid1,'%5s %17s %10s %10s %10s %10s %10s %10s %10s\n','%iobj','  centroid_xy  ','Inten_med','Area_raw','roundedness','perimeter','Equiv_Diam','Max_Diam','Solidity'); 
fprintf('%5s %17s %10s %10s %10s %10s %10s %10s %10s\n','%iobj','centroid_xy','Inten_med','Area_raw','roundedness','perimeter','Equiv_Diam','Max_Diam','Solidity'); 
lind=length(indel);

Nfilt = 0;
iobj = [];
cenxy = [];
box_xyhw = []; 
RGBobj = [];
RGBgrain_whole_selected=zeros(szI(1),szI(2),3);
Graygrain_whole_selected=zeros(szI(1),szI(2),1);
RGBom = [];
Grobj=[];
Graygrain2_scaled = [];
objN = [];
objco = [];
Igrm = [];
els2 = [];
perim = [];
sED2 = [];
sSOL2 = [];
rmaxO = [];
dmaxO = [];
permx = {};
permy = {};
lprm = [1];
subimg = {};

for idx = 1:length(indel)
    grain = false(size(bw));                        % IS THE CURRENT (idx'th) OBJECT 
    grain(cc.PixelIdxList{indel(idx)}) = true;
    contN2 = bwperim(grain,8);
    
    Lc=bwlabel(grain);
    stats = regionprops(Lc,'all');
    sA=stats(1).Area;
    scen=stats(1).Centroid;    
    sBB=stats(1).BoundingBox; % [1.5000 357.5000 8 16]
    sSAI=stats(1).SubarrayIdx; % {[358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373]  [2 3 4 5 6 7 8 9]}
    sMjA=stats(1).MajorAxisLength; 
    sMnA=stats(1).MinorAxisLength; 
    sEcen=stats(1).Eccentricity; 
    sFA=stats(1).FilledArea; 
    sED=stats(1).EquivDiameter; 
    sSOL=stats(1).Solidity; 
    sPeri=stats(1).Perimeter;
    sround=sMjA/sMnA;
       
    [yy,xx]=find(contN2==1);
    [yy1,xx1]=find(grain==1);        
    %plot(xx,yy,'.')
    
    if (length(xx)>=6)
        Araw=bwarea(grain);             % Area in pixel^2
        if (Araw>=m2snA)                       % Area filter (eliminate small objects)
            Nfilt = Nfilt + 1;
            plot(xx,yy,'.')
            plot(xx1,yy1,'.')
            ellp2=fit_ellipse(xx,yy,1);
            mjax2=ellp2.long_axis;
            mnax2=ellp2.short_axis;
            els=mjax2/mnax2;            
            text(ellp2.X0_in,ellp2.Y0_in,[num2str(indel(idx))])
            
            RGBgrain2=zeros(szI(1),szI(2),3);
            Graygrain2=zeros(szI(1),szI(2),1);
            
            for i = 1:length(xx1)
                RGBgrain_whole_selected(yy1(i),xx1(i),1) = I(yy1(i),xx1(i),1);
                RGBgrain_whole_selected(yy1(i),xx1(i),2) = I(yy1(i),xx1(i),2);
                RGBgrain_whole_selected(yy1(i),xx1(i),3) = I(yy1(i),xx1(i),3);
                RGBgrain2(yy1(i),xx1(i),1) = I(yy1(i),xx1(i),1);
                RGBgrain2(yy1(i),xx1(i),2) = I(yy1(i),xx1(i),2);
                RGBgrain2(yy1(i),xx1(i),3) = I(yy1(i),xx1(i),3);
                Graygrain_whole_selected(yy1(i),xx1(i),1) = I1(yy1(i),xx1(i),1);
                Graygrain2(yy1(i),xx1(i),1) = I1(yy1(i),xx1(i),1);
            end

%           figure
%           imshow(Graygrain2)

            subimg{Nfilt} = uint8(Graygrain2);
            
            R2=RGBgrain2(:,:,1);
            G2=RGBgrain2(:,:,2);
            B2=RGBgrain2(:,:,3);
            Gr=Graygrain2(:,:,1);
                        
            szR2=size(R2);
            R2sing=reshape(R2,(szR2(1)*szR2(2)),1);
            inzR2=find(R2sing>0);
            R2nz=R2sing(inzR2);
            
            szG2=size(G2);
            G2sing=reshape(G2,(szG2(1)*szG2(2)),1);
            inzG2=find(G2sing>0);
            G2nz=G2sing(inzG2);
            
            szB2=size(B2);
            B2sing=reshape(B2,(szB2(1)*szB2(2)),1);
            inzB2=find(B2sing>0);
            B2nz=B2sing(inzB2);
            
            szgr=size(Gr);
            Grsing=reshape(Gr,(szgr(1)*szgr(2)),1);
            inzgr=find(Grsing>0);
            Grnzgr=Grsing(inzgr);
            
            R2med=median(R2nz);
            R2mn=mean(R2nz);
                
            G2med=median(G2nz);
            G2mn=mean(G2nz);
              
            B2med=median(B2nz);
            B2mn=mean(B2nz);
                                                             
            Gr2med=median(Grnzgr);
            Gr2mn=mean(Grnzgr);
            Gr2pmicsq=log10(sum(Grnzgr))/Araw;
            
            %Gccm = (double(Gr2med).*Iscale)+MinV;
            %Gccssq = (double(Gr2pmicsq).*Iscale)+MinV;
            
            dobj0 = [];
    
            for i10 = 1:length(xx1)
                for j10 = 1:length(yy1)
                dobj0(i10,j10)=sqrt((xx1(i10)-xx1(j10))^2 + (yy1(i10)-yy1(j10))^2);
                end
            end
        
            dmax0 = max(max(dobj0));            
            rmax0 = dmax0/2;
            
            RGBobj(Nfilt,:) = [R2med G2med B2med];
            Grobj(Nfilt,:) = Gr2med;             
            iobj(Nfilt) = indel(idx);
            cenxy(Nfilt,:) = CentObj(indel(idx),:);
            box_xyhw(Nfilt,:) = BoxObj(indel(idx),:); 
            RGBobj_mean = mean(RGBobj(Nfilt,:));
            RGBom = [RGBom;RGBobj_mean];

            plot(cenxy(Nfilt,1),cenxy(Nfilt,2),'r*')%,'MarkerSize',15)
            
            %rectangle('position',[box_xyhw(Nfilt,1),box_xyhw(Nfilt,2),box_xyhw(Nfilt,3),box_xyhw(Nfilt,4)])           
            %fprintf(fid1,'%5d %8.3f %8.3f %3d %3d %3d  %3d %3d %3d %8.3f %8.3f %3d %10.2f %10.2f %10.5f %10.3f %10.3f %10.3f %10.3f\n',iobj(Nfilt),cenxy(Nfilt,:),RGBobjc(Nfilt,:),RGBobj(Nfilt,:),RGBobjc_mean,RGBobj_mean,Gr2med,Gr2medsc,Araw,Amic,Gr2pmicsq,Gr2pmicsq_sc,Gccm,Gccssq); 
            %--------------------------------------------------------------------
            % Area, Perimeter, Roundedness (major-to-minor-axis), Feret
            % Diameter stats = regionprops(BW,properties) 
            %--------------------------------------------------------------------
            
%===================================== STORE FEATURES ==========================  
%             sA=stats(1).Area;
%             scen=stats(1).Centroid;    
%             sBB=stats(1).BoundingBox; % [1.5000 357.5000 8 16]
%             sSAI=stats(1).SubarrayIdx; % {[358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373]  [2 3 4 5 6 7 8 9]}
%             sMjA=stats(1).MajorAxisLength; % 15.7292
%             sMnA=stats(1).MinorAxisLength; % 7.1282
%             sEcen=stats(1).Eccentricity; % 0.8914
%             sFA=stats(1).FilledArea; % 82
%             sED=stats(1).EquivDiameter; % 10.2179
%             sSOL=stats(1).Solidity; % 0.9111
%             sPeri=stats(1).Perimeter; % 36.6040  
%             sround=sMjA/sMnA;   
%==============================================================================

            objN = [objN; iobj(Nfilt)];
            objco = [objco; cenxy(Nfilt,:)];
            Igrm = [Igrm;Gr2med];
            AgN0 = [AgN0;Araw];
            els2 = [els2;els]; 
            perim = [perim;sPeri];
            sED2 = [sED2;sED];
            sSOL2 = [sSOL2;sSOL];
            rmaxO = [rmaxO;rmax0];
            dmaxO = [dmaxO;dmax0];
            
            %==============================================================
            %whos xx yy xx1 yy1
            permx = {permx;xx};
            permy = {permy;yy};
            lprm=[lprm;length(xx)]; 
            %==============================================================                 
            fprintf(fid1,'%5d %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',iobj(Nfilt),cenxy(Nfilt,:),Grobj(Nfilt,:),AgN0(Nfilt),els2(Nfilt),perim(Nfilt),sED2(Nfilt),dmaxO(Nfilt),sSOL2(Nfilt)); 
            fprintf('%5d %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',iobj(Nfilt),cenxy(Nfilt,:),Grobj(Nfilt,:),AgN0(Nfilt),els2(Nfilt),perim(Nfilt),sED2(Nfilt),dmaxO(Nfilt),sSOL2(Nfilt));                       
            %=============================================================================================
            % PRINT ORDER: Object_Number (iobj); Centroid (cenxy); 
            % Grayscale_Intensity_Median (Grobj);
            % Area_raw (AgN0); roundedness (els2); perimeter (perim);
            % Equivallent Diameter (sED2), Solidity (sSOL2)
            %=============================================================================================         
        end
    end
end

axis([-25 525 -25 525])
title 'Objects filtered based on shape and size'

fclose(fid1);

icr=find(isnan(RGBom)==0);
c5=corrcoef(RGBom(icr),Grobj(icr));

fprintf('%s %f\n','Correlation between mean RGB-plane & Grayscale intensities of the object median-coordinates:',c5(2,1))

%==========================================================================
figure(9)
subplot(2,2,1)
imshow(I)
title 'Raw RGB image'
subplot(2,2,2)
imshow(RGBgrain_whole)
%Iro=imrotate(RGBgrain_whole,180);
%imshow(Iro)
title 'RGB intensities dumped on all initial objects'
subplot(2,2,3)
imshow(RGBgrain_whole_selected)
%Irf=imrotate(RGBgrain_whole_selected,180);
%imshow(Irf)
title 'RGB intensities dumped on selected / filtered objects'

whos RGBgrain_whole_selected
whos permx permy

close all
figure(10)
subplot(1,2,1)
imshow(I)
title 'Raw RGB image'
subplot(1,2,2)
imshow(RGBgrain_whole_selected)
xco = objco(:,1);
yco = objco(:,2);
hold on
title 'Selected filtered objects'
mtit(num2str(filetag),'FontSize',20);

szRws=size(RGBgrain_whole_selected);

for i = 1:szRws(2)
    for j = 1:szRws(1)
        r=RGBgrain_whole_selected(j,i,1)/255;
        g=RGBgrain_whole_selected(j,i,2)/255;
        b=RGBgrain_whole_selected(j,i,3)/255;
        if (r==0 & g==0 & b==0)
        %fprintf('%d %d %d\n',r,g,b)
        else
            plot(i,j,'.','Color',[r g b])
        end
    end
end

for ii = 1:Nfilt
   req=rmaxO(ii);
   xc=cenxy(ii,1);
   yc=cenxy(ii,2);
   [XXf YYf]=open_circle(xc,yc,req);
   plot(XXf,YYf,'y.','Markersize',1.0)   
end

%plot(xco,yco,'yo','Markersize',25)
plot(xco,yco,'y.')

%plot(xint,yint,'o','Color',rgbint)

% figure(10)
% plot(xco,yco,'bo','Markersize',25)
% axis([0 Npix(2) 0 Npix(1)])

% Nfilt = length(indel);

figure(11)
subplot(1,2,1)
imshow(I1)
title 'Original Grayscale image'
%axis equal
subplot(1,2,2)
imagesc(Graygrain_whole_selected);colormap('gray')
title 'Grayscale intensities dumped on selected filtered objects'
axis equal
axis([0.5 518 0 500])
mtit(num2str(filetag),'FontSize',20);

fprintf('%d %s %d\n',Nfilt,'Objects filtered from: ',Nobj);

% figure(12)
% imagesc(Graygrain_whole_selected);colormap('gray')
% title 'Grayscale intensities dumped on selected filtered objects'
% axis equal
% axis([0.5 518 0 500])

end

%========================= END OF MAIN FUNCTION ===========================

function [XXf YYf] = open_circle(xo,yo,rad)
    XXf = [];
    YYf = [];
    deg2rad = pi/180.0;    
    for th = 0:2.5:360
        th_rad = th*deg2rad;
        xc = rad*cos(th_rad);
        yc = rad*sin(th_rad);
        XXf = [XXf;(xc+xo)];
        YYf = [YYf;(yc+yo)];
    end
end

function ellipse_t = fit_ellipse( x,y,axis_handle )
%
% fit_ellipse - finds the best fit to an ellipse for the given set of points.
%
% Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
%
% Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
%           axis_handle - optional. a handle to an axis, at which the estimated ellipse
%                         will be drawn along with it's axes
%
% Output:   ellipse_t - structure that defines the best fit to an ellipse
%                       a           - sub axis (radius) of the X axis of the non-tilt ellipse
%                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse
%                       phi         - orientation in radians of the ellipse (tilt)
%                       X0          - center at the X axis of the non-tilt ellipse
%                       Y0          - center at the Y axis of the non-tilt ellipse
%                       X0_in       - center at the X axis of the tilted ellipse
%                       Y0_in       - center at the Y axis of the tilted ellipse
%                       long_axis   - size of the long axis of the ellipse
%                       short_axis  - size of the short axis of the ellipse
%                       status      - status of detection of an ellipse
%
% Note:     if an ellipse was not detected (but a parabola or hyperbola), then
%           an empty structure is returned

% =====================================================================================
%                  Ellipse Fit using Least Squares criterion
% =====================================================================================
% We will try to fit the best ellipse to the given measurements. the mathematical
% representation of use will be the CONIC Equation of the Ellipse which is:
%
%    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
%
% The fit-estimation method of use is the Least Squares method (without any weights)
% The estimator is extracted from the following equations:
%
%    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
%
%    where:
%       A   - is the vector of parameters to be estimated (a,b,c,d,e)
%       x,y - is a single measurement
%
% We will define the cost function to be:
%
%   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
%            = (X*A+f_c)'*(X*A+f_c)
%            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
%
%   where:
%       g_c(x_c,y_c;A) - vector function of ALL the measurements
%                        each element of g_c() is g(x,y;A)
%       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
%       f_c            - is actually defined as ones(length(f),1)*f
%
% Derivation of the Cost function with respect to the vector of parameters "A" yields:
%
%   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
%
% Which yields the estimator:
%
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
%
% NOW, all that is left to do is to extract the parameters from the Conic Equation.
% We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
%
%    Recall the conic representation of an ellipse:
%
%       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
%
% We will check if the ellipse has a tilt (=orientation). The orientation is present
% if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
% tilt of the ellipse.
%
% If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
% we will remove the tilt of the ellipse so as to remain with a conic representation of an
% ellipse without a tilt, for which the math is more simple:
%
% Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
%
% We will remove the orientation using the following substitution:
%
%   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
%
%   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
%
%   where:      c = cos(phi)    ,   s = sin(phi)
%
%   and simplify...
%
%       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
%           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
%
%   The orientation is easily found by the condition of (B_new=0) which results in:
%
%   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
%
%   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
%   all the other constants A`,C`,D`,E` can be found.
%
%   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
%   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c
%   C` = A*s^2 + B*c*s + C*c^2
%
% Next, we want the representation of the non-tilted ellipse to be as:
%
%       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
%
%       where:  (X0,Y0) is the center of the ellipse
%               a,b     are the ellipse "radiuses" (or sub-axis)
%
% Using a square completion method we will define:
%
%       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
%
%       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
%                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
%
%       which yields the transformations:
%
%           X0  =   -D`/(2*A`)
%           Y0  =   -E`/(2*C`)
%           a   =   sqrt( abs( F``/A` ) )
%           b   =   sqrt( abs( F``/C` ) )
%
% And finally we can define the remaining parameters:
%
%   long_axis   = 2 * max( a,b )
%   short_axis  = 2 * min( a,b )
%   Orientation = phi
%
%

% initialize
orientation_tolerance = 1e-3;

% empty warning stack
warning( '' );

% prepare vectors, must be column vectors
x = x(:);
y = y(:);

% remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
mean_x = mean(x);
mean_y = mean(y);
x = x-mean_x;
y = y-mean_y;

% the estimation for the conic equation of the ellipse
X = [x.^2, x.*y, y.^2, x, y ];
a = sum(X)/(X'*X);

% check for warnings
if ~isempty( lastwarn )
    disp( 'stopped because of a warning regarding matrix inversion' );
    ellipse_t = [];
    return
end

% extract parameters from the conic equation
[a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5) );

% remove the orientation from the ellipse
if ( min(abs(b/a),abs(b/c)) > orientation_tolerance )
    
    orientation_rad = 1/2 * atan( b/(c-a) );
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    [a,b,c,d,e] = deal(...
        a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
        0,...
        a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
        d*cos_phi - e*sin_phi,...
        d*sin_phi + e*cos_phi );
    [mean_x,mean_y] = deal( ...
        cos_phi*mean_x - sin_phi*mean_y,...
        sin_phi*mean_x + cos_phi*mean_y );
else
    orientation_rad = 0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
end

% check if conic equation represents an ellipse
test = a*c;
switch (1)
    case (test>0),  status = '';
    case (test==0), status = 'Parabola found';  %warning( 'fit_ellipse: Did not locate an ellipse: Parabola found' );
    case (test<0),  status = 'Hyperbola found'; %warning( 'fit_ellipse: Did not locate an ellipse: Hyperbola found' );
        if (test>0)
            disp(status)
            return
        end
end

% if we found an ellipse return it's data
if (test>0)
    
    % make sure coefficients are positive as required
    if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end
    
    % final ellipse parameters
    X0          = mean_x - d/2/a;
    Y0          = mean_y - e/2/c;
    F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
    [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );
    long_axis   = 2*max(a,b);
    short_axis  = 2*min(a,b);
    
    % rotate the axes backwards to find the center point of the original TILTED ellipse
    R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
    P_in        = R * [X0;Y0];
    X0_in       = P_in(1);
    Y0_in       = P_in(2);
    
    % pack ellipse into a structure
    ellipse_t = struct( ...
        'a',a,...
        'b',b,...
        'phi',orientation_rad,...
        'X0',X0,...
        'Y0',Y0,...
        'X0_in',X0_in,...
        'Y0_in',Y0_in,...
        'long_axis',long_axis,...
        'short_axis',short_axis,...
        'status','' );
else
    % report an empty structure
    ellipse_t = struct( ...
        'a',[],...
        'b',[],...
        'phi',[],...
        'X0',[],...
        'Y0',[],...
        'X0_in',[],...
        'Y0_in',[],...
        'long_axis',[],...
        'short_axis',[],...
        'status',status );
end

% check if we need to plot an ellipse with it's axes.
if (nargin>2) & ~isempty( axis_handle ) & (test>0)
    
    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    % the axes
    ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
    horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
    hold_state = get( axis_handle,'NextPlot' );
    set( axis_handle,'NextPlot','add' );
    plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
    plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    set( axis_handle,'NextPlot',hold_state );
end
end


%MTIT		creates a major title in a figure with many axes
%
%		MTIT
%		- creates a major title above all
%		  axes in a figure
%		- preserves the stack order of
%		  the axis handles
%
%SYNTAX
%-------------------------------------------------------------------------------
%		P = MTIT(TXT,[OPT1,...,OPTn])
%		P = MTIT(FH,TXT,[OPT1,...,OPTn])
%
%INPUT
%-------------------------------------------------------------------------------
%    FH	:	a valid figure handle		[def: gcf]
%   TXT	:	title string
%
% OPT	:	argument
% -------------------------------------------
%  xoff	:	+/- displacement along X axis
%  yoff	:	+/- displacement along Y axis
%  zoff	:	+/- displacement along Z axis
%
%		title modifier pair(s)
% -------------------------------------------
%   TPx	:	TVx
%		see: get(text) for possible
%		     parameters/values
%
%OUTPUT
%-------------------------------------------------------------------------------
% par	:	parameter structure
%  .pos :	position of surrounding axis
%   .oh	:	handle of last used axis
%   .ah :	handle of invisible surrounding axis
%   .th :	handle of main title
%
%EXAMPLE
%-------------------------------------------------------------------------------
%	subplot(2,3,[1 3]);		title('PLOT 1');
%	subplot(2,3,4); 		title('PLOT 2');
%	subplot(2,3,5); 		title('PLOT 3');
%	axes('units','inches',...
%	     'color',[0 1 .5],...
%	     'position',[.5 .5 2 2]);	title('PLOT 41');
%	axes('units','inches',...
%	     'color',[0 .5 1],...
%	     'position',[3.5 .5 2 2]);	title('PLOT 42');
%	shg;
%	p=mtit('the BIG title',...
%	     'fontsize',14,'color',[1 0 0],...
%	     'xoff',-.1,'yoff',.025);
% % refine title using its handle <p.th>
%	set(p.th,'edgecolor',.5*[1 1 1]);

% created:
%	us	24-Feb-2003		/ R13
% modified:
%	us	24-Feb-2003		/ CSSM
%	us	06-Apr-2003		/ TMW
%	us	13-Nov-2009 17:38:17

%-------------------------------------------------------------------------------
function par=mtit(varargin)

defunit='normalized';
if	nargout
    par=[];
end

% check input
if	nargin < 1
    help(mfilename);
    return;
end
if	isempty(get(0,'currentfigure'))
    disp('MTIT> no figure');
    return;
end

vl=true(size(varargin));
if	ischar(varargin{1})
    vl(1)=false;
    figh=gcf;
    txt=varargin{1};
elseif	any(ishandle(varargin{1}(:)))		&&...
        ischar(varargin{2})
    vl(1:2)=false;
    figh=varargin{1};
    txt=varargin{2};
else
    error('MTIT> invalid input');
end
vin=varargin(vl);
[off,vout]=get_off(vin{:});

% find surrounding box
ah=findall(figh,'type','axes');
if	isempty(ah)
    disp('MTIT> no axis');
    return;
end
oah=ah(1);

ou=get(ah,'units');
set(ah,'units',defunit);
ap=get(ah,'position');
if	iscell(ap)
    ap=cell2mat(get(ah,'position'));
end
ap=[	min(ap(:,1)),max(ap(:,1)+ap(:,3)),...
    min(ap(:,2)),max(ap(:,2)+ap(:,4))];
ap=[	ap(1),ap(3),...
    ap(2)-ap(1),ap(4)-ap(3)];

% create axis...
xh=axes('position',ap);
% ...and title
th=title(txt,vout{:});
tp=get(th,'position');
set(th,'position',tp+off);
set(xh,'visible','off','hittest','on');
set(th,'visible','on');

% reset original units
ix=find(~strcmpi(ou,defunit));
if	~isempty(ix)
    for	i=ix(:).'
        set(ah(i),'units',ou{i});
    end
end

% ...and axis' order
uistack(xh,'bottom');
axes(oah);				%#ok

if	nargout
    par.pos=ap;
    par.oh=oah;
    par.ah=xh;
    par.th=th;
end
end
%-------------------------------------------------------------------------------
function	[off,vout]=get_off(varargin)

% search for pairs <.off>/<value>

off=zeros(1,3);
io=0;
for	mode={'xoff','yoff','zoff'};
    ix=strcmpi(varargin,mode);
    if	any(ix)
        io=io+1;
        yx=find(ix);
        ix(yx+1)=1;
        off(1,io)=varargin{yx(end)+1};
        varargin=varargin(xor(ix,1));
    end
end
vout=varargin;
end
%-------------------------------------------------------------------------------

function [minVmv,maxVmv]=returnMinMax(finp,datmm);

%datmm='Min_Max_Table_example_phase_data.format';

minmax=importdata(datmm);
minmaxV=minmax.data;
rows=minmax.textdata;
header=rows(1,:);
% dirname=rows(2:length(rows),1);
% fname=rows(2:length(rows),2);

fname = rows(2:length(rows),1);

minV=minmaxV(:,1);
maxV=minmaxV(:,2);

%whos dirname fname minV maxV

% dinp=input('Enter dirname:\n','s');
% finp=input('Enter filename:\n','s');

flag=0;
for i = 1:length(fname)
    %dirstr=dirname{i};
    fnamestr=fname{i};
    %whos str1 dinp
    %match1=strcmp(dirstr,dinp);
    match=strcmp(fnamestr,finp);
    disp(match);
    flag=0;
    if (match == 1) % & match2 == 1)
        minVmv=minV(i);
        maxVmv=maxV(i);
        %fprintf('%20s  %5d  %15s  %15s  %10.3f  %10.3f\n','Match found at i=',i,dirstr,fnamestr,minVmv,maxVmv)
        flag=1;
        break;
    end
end

    if (flag==0)
    fprintf('Data not found for %20s\n',fnamestr)
    return
    end
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
    newimA{i}=z(30:429,30:429);      % Cropping done here (user have to reset appropriately)
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

