function [dp, pvalue]=getDprime(data1,data2)
if ~isvector(data1) || ~isvector(data2)
    error('The data should be vectors');
end
if size(data1,1)~=1
    data1=transpose(data1);
end
if size(data2,1)~=1
    data2=transpose(data2);
end
mu1=mean(data1);
mu2=mean(data2);
ee1=std(data1);
ee2=std(data2);

dp=(mu1-mu2)/sqrt((ee1^2+ee2^2)/2);
if nargout==2
    data=[data1 data2];
    nn1=length(data1);
    nn2=length(data2);
    Ntrl=200;
    dp_new=zeros(1, Ntrl);
    for ii=1:Ntrl
        idxNew=randi(nn1+nn2, 2, nn1);
        dataNew1=data(idxNew(1,:));
        dataNew2=data(idxNew(2,:));
        mu1=mean(dataNew1);
        mu2=mean(dataNew2);
        ee1=std(dataNew1);
        ee2=std(dataNew2);
        dp_new(ii)=(mu1-mu2)/sqrt((ee1^2+ee2^2)/2);    
    end
    
    pvalue=sum(dp_new>=dp)/Ntrl;
    if pvalue>0.5
        pvalue=1-pvalue;
    end
end
