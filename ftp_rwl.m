function ftp_rwl(path,file)

    full_path=char(path);
    filename=char(file);
    website=full_path(7:23);
    directory=full_path(24:end);
    
    ftpobj=ftp(website);
    cd(ftpobj,directory);
    mget(ftpobj,filename,'./rwls_raw');

end
