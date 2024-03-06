
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
            for n2=n1+1:n-2
                for m2=n2+2:n
                    v2=v(n2:m2);
                    i2=i(n2:m2);
                    X=zeros(3,3);
                    X(1,1)=length(v2);
                    X(1,2)=sum(v2);
                    X(1,3)=sum(i2);
                    X(2,1)=X(1,2);
                    X(2,2)=sum(v2.*v2);
                    X(2,3)=sum(v2.*i2);
                    X(3,1)=X(1,3);
                    X(3,2)=X(2,3);
                    X(3,3)=sum(i2.*i2);
                    Y=zeros(3,1);
                    term1=log(K-(F*v2)-i2);
                    Y(1,1)=sum(term1);
                    Y(2,1)=sum(v2.*term1);
                    Y(3,1)=sum(i2.*term1);
                    BCD=X\Y;
                    BCD=exp(BCD);
                    B=BCD(1);
                    C=BCD(2);
                    D=BCD(3);
                        for n3=n1+1:n-1
                            for m3=n3+1:n
                                v3=v(n3:m3);
                                i3=i(n3:m3);
                                term1=(v3*log(C))+(i3*log(D));
                                term2=log((K-(F*v3))-(B*(C.^v3).*(D.^i3))-(i3));
                                X=zeros(2,2);
                                X(1,1)=length(v3);
                                X(1,2)=sum(term1);
                                X(2,1)=X(1,2);
                                X(2,2)=sum(term1.*term1);
                                Y=zeros(2,1);
                                Y(1,1)=sum(term2);
                                Y(2,1)=sum(term2.*term1);
                                Eg=X\Y;
                               E=exp(Eg(1));
                               g=Eg(2);
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