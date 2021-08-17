clear
close all

flist=dir('*.rwl');
baddogs={};
k=1;
for i =1:length(flist);
    filename=flist(i).name;
    disp('------------------------------------')
    disp(['now testing ' filename ])
    disp('------------------------------------')
    pause(.1);
    
    try 
        [X,yrX,nms,T]=rwl2tsm(['./' filename]);        
    catch
        
        disp([filename ' will not work'])
        baddogs{k}=filename;
        k=k+1;
        
    end
end
%%
disp('------------------------------------')
disp(['End of list - SUCCESS!!!' ])
disp('------------------------------------')
disp(['The list below are problem files:' ])
disp(char(baddogs));


