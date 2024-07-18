function [chromo,temp1,remain_pop] = elitism( pop,combine_chromo2,f_num,x_num )
[pop2,temp]=size(combine_chromo2);
chromo_rank=zeros(pop2,temp);
chromo=[];
[~,index_rank]=sort(combine_chromo2(:,f_num+x_num+1));
for i=1:pop2
    chromo_rank(i,:)=combine_chromo2(index_rank(i),:);
end
max_rank=max(combine_chromo2(:,f_num+x_num+1));
prev_index=0;
for i=1:max_rank
    current_index=max(find(chromo_rank(:,(f_num+x_num+1))==i));
    if(current_index>pop)
        remain_pop=pop-prev_index;
        temp1=chromo_rank(((prev_index+1):current_index),:);
        return;
    elseif(current_index<=pop)
        chromo(((prev_index+1):current_index),:)=chromo_rank(((prev_index+1):current_index),:);
    end
    prev_index = current_index;
end
end
