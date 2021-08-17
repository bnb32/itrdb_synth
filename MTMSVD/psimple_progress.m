function psimple_progress(cdir,i,N);
    w=50;
    if i==0;
        delete([cdir,'/count*.txt']);
    else
        fname=[cdir,'/count',num2str(i),'.txt'];
        f=fopen(fname,'a');
        fclose(f);
        a=dir([cdir,'/count*.txt']);
        out=size(a,1);
        percent=out/N*100;
        progress=sprintf('%3.0f%',percent);
        %fprintf([progress,'%%\n']);
        fprintf([repmat(char(8),1,w+9),progress,'%%[',repmat('=',1,round(percent*w/100)), '>', repmat(' ',1,w-round(percent*w/100)),']\n']);
    end
end    
