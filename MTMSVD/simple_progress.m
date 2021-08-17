function simple_progress(idx,num,div);
    for i=1:div
        bins(i)=ceil(i*num/div);
    end
    if any(idx==bins);
        disp([num2str(idx) '/' num2str(num)]);
    end	
end    
