
clear all
clc
input2;
x=load('two_diode_valid_results.txt');

  [m, n]=size(x);
  kk=1;
  rmse_min=zeros(1,m);
  
   for i=1:m
      Iph=x(i,1);
      I01=x(i,2);
      I02=x(i,3);
      a1=x(i,4);
      a2=x(i,5);
      Rs=x(i,6);
      Rp=x(i,7);
      sim('pvmodelling.slx');
      xx=load('sree.mat');
      i1=xx.ii';
      i1=i1(:,2);
      data=load('1.txt');
      i2=data(:,2);
                 
      rmse=sqrt(sum((i1-i2).*(i1-i2))/length(i1));
      rmse_min(kk)=rmse;
      kk=kk+1;
      current=[i1   i2]
  end
   
    rmse_min
    
   