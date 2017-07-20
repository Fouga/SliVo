function num = find_equi_number(total_num)

    if total_num < 1500 
        num = total_num-1;
    elseif total_num > 1500 && total_num < 5000
        num = 1500;
    elseif total_num > 5000 && total_num < 10000
        num = floor(total_num*1/3);         
    elseif total_num > 10000 && total_num < 40000
        num = floor((total_num)*1/4);
    elseif total_num > 40000
        num = 10000;
    end



%     if total_num < 1000 
%         num = total_num;
%     elseif total_num > 1000 || total_num < 3000
%         num = floor(total_num*1/2);
%     elseif total_num > 3000 || total_num < 6000
%         num = floor(total_num*1/3);  
%     elseif total_num > 6000 || total_num < 10000
%         num = floor(total_num*1/4);         
%     elseif total_num > 10000 || total_num < 40000
%         num = floor((total_num).^(3/4));
%     elseif total_num > 40000
%         num = 10000;
%     end