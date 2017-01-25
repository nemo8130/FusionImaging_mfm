
function [I1sc]=scaleIntenMM(file)

%file='aptes_010_obj.tif.obj';

datmm='Min_Max_Table_example_phase_data.format.2';

[MinV,MaxV] = returnMinMax(file,datmm);
fprintf('%s %f %s %f\n','MinV:',MinV,'MaxV:',MaxV)

dtall=importdata(file);
dt=dtall.data;
I1=dt(:,4);

Iscale = ((MaxV-MinV)/255);

fprintf('%s %s %d %s %d\n','Before Scaling:','MinPx:',min(min(I1)),'MaxPx:',max(max(I1)))
I1sc = (double(I1).*Iscale)+MinV;
fprintf('%s %s %f %s %f\n','After Scaling:','MinScPx:',min(min(I1sc)),'MaxScPx:',max(max(I1sc)))

subplot(2,1,1);plot(I1);ylabel 'Intensity';xlabel 'coordinates';subplot(2,1,2);plot(I1sc,'r');ylabel 'Scaled Intensities';xlabel 'coordinates';

outfile=strcat(file,'.scaled');
fid1=fopen(outfile,'w');
fprintf(fid1,'%5s %17s %10s %10s %10s %10s %10s %10s %10s\n','%iobj','  centroid_xy  ','Scl_Inten','Area_raw','roundedness','perimeter','Equiv_Diam','Max_Diam','Solidity'); 

    for i = 1:length(dt)
        fprintf(fid1,'%5d %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',dt(i,1),dt(i,2),dt(i,3),I1sc(i),dt(i,5),dt(i,6),dt(i,7),dt(i,8),dt(i,9),dt(i,10)); 
    end

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

farr = strsplit(finp,'.');
finpn=farr(1);

flag=0;
for i = 1:length(fname)
    %dirstr=dirname{i};
    fnamestr=fname{i};
    fnarr = strsplit(fnamestr,'.');
    fnloop=fnarr(1);
    %whos str1 dinp
    %match1=strcmp(dirstr,dinp);
    match=strcmp(fnloop,finpn);
    disp(match);%disp(match2)
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

