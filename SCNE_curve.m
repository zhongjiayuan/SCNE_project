
clear;
clc;
close all;
[data_entropy,empty] = xlsread('KIRC_entropy_matrix.xlsx');

pnum=[198,41,64,50];
data_entropy_size=size(data_entropy);
data_entropy(data_entropy>10^3)=mean(mean(data_entropy));

n=1;
for t=1:length(pnum)
    count=0;
    for c=n:sum(pnum(1:t))   %psize(3)
        count=count+1;
        tmp_entropy=sort(data_entropy(:,c),'descend');
        aver_entropy(t,count)=mean(tmp_entropy(1:floor(data_entropy_size(1)*0.05)));
    end
    case_result(t)=mean(aver_entropy(t,1:pnum(t)));
    n=1+sum(pnum(1:t));
end

plot(1:4,case_result,'r-*','LineWidth',3);








