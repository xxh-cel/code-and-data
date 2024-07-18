clear all
clc
tic;
ea=1;
U=50;
DD=1;
ea=1;
HV=[];
xx=[];
ddd=[];
zt=[];
kss=0;
seta=1;
zhibiao=zeros(30,2);
fun='DTLZ7';
pop=66;
gen=50;
maxreal=150;
GZMIN=[];
QW=1;
seta=0;
zhuangtai=1;
f_num=3;
x_num=6;
 save M.mat  f_num
x_min=zeros(1,x_num);
x_max=[150,150,150,150,150,150];
pc=1;
pm=1/6;
yita1=30;
yita2=20;
[Z,KK] = UniformPoint(pop,f_num);
pdate = cell(1, 20);
for i = 1:80 
   a = [
        floor(rand() * 33) + 1,
        floor(rand() * 33) + 34,
        floor(rand() * 53) + 67
    ];
    pdate{i} = a;
end
pTable = cell2table(pdate);
writetable(pTable, '3mubiaobugudingriqi.csv');
ig=1;
for ff =1:30
    save ff.mat  ff
    dat=pdate(ff)
    save p.mat  dat
    chromo=floor(initialize(pop,f_num,x_num,x_min,x_max,fun));
    save chromo.mat chromo
    Rpath = 'E:\R-4.3.1\bin\Rscript.exe'; 
    RscriptFileName = 'F:\dssa\code\mr.R'; 
    RunRcode(RscriptFileName, Rpath)
    dd=load('FUNV.mat');
    chromo(:,x_num+1:x_num+f_num)=dd.FunctionValue;
    spo=sum(chromo(:,x_num+1:x_num+f_num),1);
    spo=spo/pop;
    [F1,chromo_non]=non_domination_sort(pop,chromo,f_num,x_num);
    i=1;
    NoZ = size(Z,1);
    FunctionValue=chromo(:,x_num+1:x_num+f_num);
    [N,M] = size(FunctionValue);
    Zmin = min(FunctionValue,[],1);	
    Extreme = zeros(1,M);          
    w = zeros(M)+0.000001+eye(M);   
    for i = 1 : M                  
        [~,Extreme(i)] = min(max(FunctionValue./repmat(w(i,:),N,1),[],2));
    end
    Hyperplane = FunctionValue(Extreme,:)\ones(M,1);	
    a = 1./Hyperplane;             
    if any(isnan(a))
        a = max(FunctionValue,[],1)';
    end
    FunctionValue = (FunctionValue-repmat(Zmin,N,1))./(repmat(a',N,1)-repmat(Zmin,N,1));	
    Distance = zeros(N,NoZ);      
    normZ = sum(Z.^2,2).^0.5;      
    normF = sum(FunctionValue.^2,2).^0.5;   
    for i = 1 : N
        for j = 1 : NoZ  
           S1 = normF(i);
           S2 = normZ(j);
           cos(i,j)= sum(FunctionValue(i,:).*Z(j,:),2)./S2./S1;
        end
     end
     d11=zeros(N,NoZ);
     D1=zeros(N,NoZ);
     for i=1:N
         for j = 1 : NoZ
            S1 = normF(i);
            d11(i,j)=S1*cos(i,j);
            d2(i,j)=S1*sqrt(1-cos(i,j).^2);
            D1(i,j)=d11(i,j)+seta*d2(i,j);  
        end
     end
     [d0,pi]=min(d2',[],1);
     for i= 1:N
        d0(i)=D1(i,pi(i));
     end
     rho = zeros(1,NoZ);
     for i = 1 : N
        rho(pi(i)) = rho(pi(i))+1;
     end
     chromo_offspring=feicross_mutation1(chromo,f_num,x_num,x_min,x_max,pc,pm,yita1,yita2,fun ,chromo_non,d0,pi,rho,pop);
     chromo_parent=chromo;
     chromo_offspring=chromo_offspring(:,1:x_num);
     save chromo_offspring.mat chromo_offspring
     Rpath = 'E:\R-4.3.1\bin\Rscript.exe'; 
     RscriptFileName1 = 'F:\dssa\code\mr1.R'; 
     RunRcode(RscriptFileName1, Rpath)
     dd1=load('offFUNV.mat');
     chromo_offspring(:,x_num+1:x_num+f_num)=dd1.FunctionValue ; 
     while ea<gen
   save gen.mat  ea
    [pop_parent,~]=size(chromo);
    [pop_offspring,~]=size(chromo_offspring);
    combine_chromo(1:pop_parent,1:(f_num+x_num))=chromo(:,1:(f_num+x_num));
    combine_chromo((pop_parent+1):(pop_parent+pop_offspring),1:(f_num+x_num))=chromo_offspring(:,1:(f_num+x_num));
    [pop2,~]=size(combine_chromo);
    [F2,combine_chromo1]=non_domination_sort(pop2,combine_chromo,f_num,x_num);
    [temp,temp1,ks]=elitism(pop,combine_chromo1,f_num,x_num);  
    temp_size=size(temp,1);
    if temp_size==pop
       chromo=temp;
       if(kss==0)
         seta=5;
         kss=1;
       else
           seta=seta;
        end
    else
       if size(temp)==0
         [Choose,d,pi,rho,seta,xx] = F_choose11(temp,temp1(:,x_num+1:x_num+f_num),ks,Z,zhuangtai,ea,DD,seta,gen,xx);
         temp1=temp1(Choose,:);
         chromo=[temp;temp1];
       else
          [Choose,d,pi,rho,seta,xx] = F_choose11(temp(:,x_num+1:x_num+f_num),temp1(:,x_num+1:x_num+f_num),ks,Z,zhuangtai,ea,DD,seta,gen,xx);
          temp1=temp1(Choose,:);
          chromo=[temp;temp1];
       end
    end
    spo11=sum(chromo(:,x_num+1:x_num+f_num),1);
    pjl=chromo(:,x_num+1:x_num+f_num);
    spo1=sum(pjl(:,1:f_num),1);
    spo1=spo1/pop;
    zhuangtai=sum((spo-spo1).^2,2).^0.5/sum(spo.^2,2).^0.5;
    zt(end+1)=zhuangtai
    spo=spo1; 
    NoZ = size(Z,1); 
    FunctionValue=chromo(:,x_num+1:x_num+f_num);
    [N,M] = size(FunctionValue);  
    length(chromo)
    Zmin = min(FunctionValue,[],1);
    Extreme = zeros(1,M);          
    w = zeros(M)+0.000001+eye(M);  
    for i = 1 : M                 
        [~,Extreme(i)] = min(max(FunctionValue./repmat(w(i,:),N,1),[],2));
    end
    Hyperplane = FunctionValue(Extreme,:)\ones(M,1);	
    a = 1./Hyperplane;           
   
        if any(isnan(a))
        a = max(FunctionValue,[],1)';
        end
    a(a < 0.001* max(FunctionValue,[],1)') = max(FunctionValue(:,a < 0.001* max(FunctionValue,[],1)'),[],1); 
    FunctionValue = (FunctionValue-repmat(Zmin,N,1))./(repmat(a',N,1)-repmat(Zmin,N,1));	
    Distance = zeros(N,NoZ);       
    normZ = sum(Z.^2,2).^0.5;       
    normF = sum(FunctionValue.^2,2).^0.5;   
         for i = 1 : N
            for j = 1 : NoZ
               S1 = normF(i);
               S2 = normZ(j);
               cos(i,j)= sum(FunctionValue(i,:).*Z(j,:),2)./S2./S1;
            end
          end
    d11=zeros(N,NoZ);
     D1=zeros(N,NoZ);
    for i=1:N
         for j = 1 : NoZ
              S1 = normF(i);
             d11(i,j)=S1*cos(i,j);
             d2(i,j)=S1*sqrt(1-cos(i,j).^2);
             D1(i,j)=d11(i,j)+seta*d2(i,j);    
        end
    end
    [d,pi]=min(d2',[],1);
    for i= 1:N
        d(i)=D1(i,pi(i));
    end
    DD=abs(sum(d)/pop-sum(d0)/pop)/(sum(d0)/pop)
    ddd(end+1)=DD
    d0=d;
    rho = zeros(1,NoZ);
    for i = 1 : N
        rho(pi(i)) = rho(pi(i))+1;
    end
    [F1,chromo_non]=non_domination_sort(pop,chromo,f_num,x_num);
    chromo_offspring=feicross_mutation1(chromo,f_num,x_num,x_min,x_max,pc,pm,yita1,yita2,fun ,chromo_non,d,pi,rho,pop);
    chromo_offspring=chromo_offspring(:,1:x_num);
    save chromo_offspring.mat chromo_offspring
    Rpath = 'E:\R-4.3.1\bin\Rscript.exe'; 
    RscriptFileName1 = 'F:\dssa\code\mr1.R'; 
    RunRcode(RscriptFileName1, Rpath)
    dd1=load('offFUNV.mat');
    chromo_offspring(:,x_num+1:x_num+f_num)=dd1.FunctionValue ; 
   [F1,chromo_non]=non_domination_sort(pop,chromo_offspring,f_num,x_num);
   ea=ea+1;
    end
    end






