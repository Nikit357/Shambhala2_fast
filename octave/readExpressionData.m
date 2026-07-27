function data=readExpressionData(filename,isLog)
data=struct('GeneList',{},'SamplesName',{},'Samples',{});
infile=fopen(filename,'r');
headLine=fgets(infile);
headName=textscan(headLine,'%s');
N=size(headName{1},1);
strDataFormat='%s';
Cells=0;
iSamples=0;
for i=2:N
    iSamples=iSamples+1;
    Cells(iSamples)=i;
    strDataFormat = [strDataFormat, '%f64'];
end
inData=textscan(infile,strDataFormat);
fclose(infile);
volSample=size(inData{1},1);
data(1).Samples=zeros(volSample,N-1);
data.GeneList=inData{1};
for i=1:N-1
    if(strcmp(isLog,'log2'))
        data.Samples(:,i)=log2(inData{1,Cells(i)}+1);
    else
        data.Samples(:,i)=inData{1,Cells(i)};
    end
    data.SamplesName{i}=headName{1}{Cells(i)};
end
    