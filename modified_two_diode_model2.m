
clear all
clc
input2;
data=load('1.txt');
v=data(:,1);
i=data(:,2);
n=size(v);
iter=1;
r1=zeros(6798,7);
r2=zeros(6798,6);
for n1=1:n-2
    for m1=n1+1:n-1
        v1=v(n1:m1);
        i1=i(n1:m1);
        X=zeros(2,2);
        X(1,1)=length(v1);
        X(1,2)=-sum(v1);
        X(2,1)=-X(1,2);
        X(2,2)=-sum(v1.*v1);
        Y=zeros(2,1);
        Y(1,1)=sum(i1);
        Y(2,1)=sum(v1.*i1);
        KF=X\Y;
        K=KF(1);
        F=KF(2);
            

        %% Add code like above i.e., B, C & D parameter extraction

        %% Add code like above i.e., E & g parameter extraction

                                         A=K-B-E;
                                         zeta1=1/(log(C)*Vt);
                                         zeta2=(zeta1)/g;
                                         rs=zeta1*Vt*log(D);
                                         rp=(1/F)-rs;
                                         Iph=(A/rp)*(rp+rs);
                                         I01=(B/rp)*(rp+rs);
                                         I02=(E/rp)*(rp+rs);
                                parame=[Iph I01 I02 zeta1 zeta2 rs rp];
                                
                                if(min(parame)>0 && isreal(parame)==1)
                                    if(zeta1>=1&& zeta2>=1)
                                         disp(parame);
                                            r1(iter,:)=parame;
                                            r2(iter,:)=[n1 m1 n2 m2 n3 m3];
                                            iter=iter+1;
                                    end
                                end

                                end
                                
                             end
                        end
                end
            end
end
    r3=[r1   r2];
save('two_diode_valid_results.txt','r3','-ascii');
disp(iter);