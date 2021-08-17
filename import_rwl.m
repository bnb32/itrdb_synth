function [results,imp_fail,badrwl,rwlerr]=import_rwl(local_path,filename);
    results='';
    imp_fail=false;
    badrwl='';
    rwlerr='';
    try 
	results=rwl2tsm([local_path filename]);
    catch
        %error('STOP!')
        pause(.5)
        badrwl=filename;
	rwlerr=lasterr;
        %badrwl=filename;
	imp_fail=true;
    end
end

