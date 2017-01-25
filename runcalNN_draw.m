% Draw corresponding objects from a pairwise Nearest-Neighbor Algo and
% replicate it for all images
% listfile='aptes_newMM.list';
% USAGE: runcalNN_draw('aptes_newMM.list')
% 
function runcalNN_draw(listfile)

load aptes_refdat.mat;

list=importdata(listfile);
file1=list(1,:);
f1=importdata(char(file1));
iitarg = [1:length(f1.data)];
xorg=f1.data(:,2);
yorg=f1.data(:,3);

icnt = zeros(length(iitarg),1);
iimap = [];

figure
hold on

for kk = 2:length(list)-1
        file2=list(kk,:);
        ff1 = char(file1);
        ff2 = char(file2);    
        fprintf('%s %s\n',ff1,ff2)
        [ii,tco,nco]=calNN_draw(ff1,ff2,0,1);
        
        for j = 1:length(iitarg)
            for k = 1:length(ii)
                if (iitarg(j)==ii(k))
                    icnt(j)=icnt(j)+1;
                end
            end
        end
        %iimap = [iimap ii];
end

title ('APTES: Corresponding Objects Identified from Nearest Neighbor Algorithm','FontSize',15)

figure
subplot(1,2,1)
imagesc(Graygrain_whole_selected)
colormap('gray'); axis equal; axis([0.5 518 0 500])
hold on
ifound=find(icnt>=(length(list)-2));
iitargf = iitarg(ifound);
fprintf('%s %d %s %d\n','Number of objects appeared in all but one snapshots',length(ifound),' out of ',length(iitarg));

%text(xorg(ifound)+1,yorg(ifound)+1,num2str(iitargf),'FontSize',15)

for kk = 2:length(list)-1
        file2=list(kk,:);
        ff1 = char(file1);
        ff2 = char(file2);    
        fprintf('%s %s\n',ff1,ff2)
        [ii,tco,nco]=calNN_draw(ff1,ff2,0,0);
        if (kk==2)
            plot(tco(:,1),tco(:,2),'yo',...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[.49 1 .63],...
                    'MarkerSize',7.5)
%             [XXf YYf]=open_circle(xm,ym,rad);
%             plot(XXf,YYf,'y.','Markersize',1.0)
%             text(xm+2,ym-5,num2str(ik),'FontSize',15,'Color','y')
        end
        for j = 1:length(iitargf)
            for k = 1:length(ii)
                if (iitargf(j)==ii(k))
                    plot(nco(k,1),nco(k,2),'bo',...
                    'MarkerEdgeColor','c',...
                    'MarkerFaceColor','m',...
                    'MarkerSize',6.0)                  
                    XX = [tco(k,1);nco(k,1)];
                    YY = [tco(k,2);nco(k,2)];
                    plot(XX,YY,'g-','LineWidth',1.5)
                end
            end
        end
end
plot(xorg(ifound),yorg(ifound),'k*')
Xbox=[150;(150+100);(150+100);150;150]
Ybox=[425;425;(425+50);(425+50);425]
plot(Xbox,Ybox,'Linewidth',1.5,'Color','w')
title ('APTES: Objects that are present in 6 out of 7 snapshots','FontSize',12)
%legend('Object centroids (All)','Corresponding objects in other snapshots','Connecting corresponding objects by the NN-algo','Location','Best')

subplot(1,2,2)
imagesc(Graygrain_whole_selected)
colormap('gray'); axis equal; axis([0.5 518 0 500])
hold on
ifound=find(icnt>=(length(list)-2));
iitargf = iitarg(ifound);
fprintf('%s %d %s %d\n','Number of objects appeared in all but one snapshots',length(ifound),' out of ',length(iitarg));

%text(xorg(ifound)+1,yorg(ifound)+1,num2str(iitargf),'FontSize',15)

for kk = 2:length(list)-1
        file2=list(kk,:);
        ff1 = char(file1);
        ff2 = char(file2);    
        fprintf('%s %s\n',ff1,ff2)
        [ii,tco,nco]=calNN_draw(ff1,ff2,0,0);
        if (kk==2)
            plot(tco(:,1),tco(:,2),'yo',...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[.49 1 .63],...
                    'MarkerSize',7.5)
%             [XXf YYf]=open_circle(xm,ym,rad);
%             plot(XXf,YYf,'y.','Markersize',1.0)
%             text(xm+2,ym-5,num2str(ik),'FontSize',15,'Color','y')
        end
        for j = 1:length(iitargf)
            for k = 1:length(ii)
                if (iitargf(j)==ii(k))
                    plot(nco(k,1),nco(k,2),'bo',...
                    'MarkerEdgeColor','c',...
                    'MarkerFaceColor','m',...
                    'MarkerSize',7.5)
                    XX = [tco(k,1);nco(k,1)];
                    YY = [tco(k,2);nco(k,2)];
                    plot(XX,YY,'g-','LineWidth',1.5)
                end
            end
        end
end
plot(xorg(ifound),yorg(ifound),'k*')
title ('A single object (boxed in the left panel) magnified')
mtit ('Identification of Corresponding Objects by the Nearest Neighbor Algorithm','FontSize',12)
legend('Object centroids (All)','Corresponding objects in other snapshots','Connecting corresponding objects by the NN-algo','Location','Best')

axis([150 250 425 475])

end
