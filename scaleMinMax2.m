
list = importdata('aptes_new.list');
coll = [];

for i = 1:length(list)
    datw=importdata(list{i});
    datn=datw.data;
    coll=[coll;datn(:,4)];
end

minI=min(coll);
maxI=max(coll);
Isccol = [];

for i = 1:length(list)
    datw=importdata(list{i});
    datn=datw.data;
    I=datn(:,4);
    Isc=(I-minI)/(maxI-minI);
    Isccol = [Isccol;Isc];
    length(Isc)
    newfile=strcat(list{i},'_minmax');
    fid1=fopen(newfile,'w');
    fprintf(fid1,'%5s %17s %10s %10s %10s %10s %10s %10s %10s\n','%iobj','  centroid_xy  ','Inten_med','Area_raw','roundedness','perimeter','Equiv_Diam','Max_Diam','Solidity');
    for j = 1:length(Isc)
        fprintf(fid1,'%5d %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',datn(j,1),datn(j,2),datn(j,3),Isc(j),datn(j,5),datn(j,6),datn(j,7),datn(j,8),datn(j,9),datn(j,10));
    end
    fclose (fid1);
end

max(Isccol)
min(Isccol)





