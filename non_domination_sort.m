function [F,chromo] = non_domination_sort( pop,chromo,f_num,x_num )
pareto_rank=1;
F(pareto_rank).ss=[];
p=[];
for i=1:pop
    p(i).n=0;
    p(i).s=[];
    for j=1:pop
        less=0;
        equal=0;
        greater=0;
        for k=1:f_num
            if(chromo(i,x_num+k)<chromo(j,x_num+k))
                less=less+1;
            elseif(chromo(i,x_num+k)==chromo(j,x_num+k))
                equal=equal+1;
            else
                greater=greater+1;
            end
        end
        if(less==0 && equal~=f_num)
            p(i).n=p(i).n+1;
        elseif(greater==0 && equal~=f_num)
            p(i).s=[p(i).s j];
        end
    end
    if(p(i).n==0)
        chromo(i,f_num+x_num+1)=1;
        F(pareto_rank).ss=[F(pareto_rank).ss i];
    end
end
while ~isempty(F(pareto_rank).ss)
    temp=[];
    for i=1:length(F(pareto_rank).ss)
        if ~isempty(p(F(pareto_rank).ss(i)).s)
            for j=1:length(p(F(pareto_rank).ss(i)).s)
                p(p(F(pareto_rank).ss(i)).s(j)).n=p(p(F(pareto_rank).ss(i)).s(j)).n - 1;
                if p(p(F(pareto_rank).ss(i)).s(j)).n==0
                    chromo(p(F(pareto_rank).ss(i)).s(j),f_num+x_num+1)=pareto_rank+1;
                    temp=[temp p(F(pareto_rank).ss(i)).s(j)];
                end
            end
        end
    end
    pareto_rank=pareto_rank+1;
    F(pareto_rank).ss=temp;
end
