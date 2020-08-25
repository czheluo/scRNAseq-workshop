%PCA FOR scRNA
x=Y_PCA(:,1:3);
[n,k]=size(x);
 x1=x;
OF=[];tm=[];
%g=input(' g=');
g=5;
k1ssw=rand(1,g)/1000;
% initial Clustering
%inits=input('Initial clustering, Evenly partitioned/=1 or most unevenly one/=2 ')
figure(1),clf
for init=1:1
    ta1=cputime;
    disp(['group, init:',' ',num2str(g),',  ',num2str(init)])
    y1=nan(n,k,g);n1=zeros(1,k);
    if init<=1
        ni=floor(n/g);n2=1;
        for i=1:g-1
            n1(i)=ni;
            y1(n2:n2+n1(i)-1,:,i)=x(n2:n2+n1(i)-1,:);
            n2=n2+n1(i);
        end
        n1(g)=n-n2+1;
        y1(n2:n,:,g)=x(n2:n,:);
    elseif init==2
        if yn==0
            x=x1;
        else
            x=x1./(ones(n,1)*std(x1));
        end
        n2=1;
        for i=1:g-1
            n1(i)=1;
            y1(n2:n2,:,i)=x(n2:n2,:);
            n2=n2+1;
        end
        n1(g)=n-n2+1;
        y1(n2:n,:,g)=x(n2:n,:);
    else
        if yn==0
            x=x1;
        else
            x=x1./(ones(n,1)*std(x1));
        end
        while min(n1)==0
            y1=nan(n,k,g);n1=zeros(1,g);gi=random('unid',g,n,1);
            for l=1:n
                y1(l,:,gi(l))=x(l,:);
                n1(gi(l))=n1(gi(l))+1;
            end
        end
    end
    r3=range(x1(:,k));
    ym=y1;
    figure(1),clf
    hold off
    cagraph
    pause(1.001)
    % calculate SSw
    SSw=zeros(k);ti=zeros(g,k);
    for i=1:g
        ss=zeros(k);
        for l=1:n
            if ~isnan(y1(l,:,i))
                ss=ss+y1(l,:,i)'*y1(l,:,i);
                ti(i,:)=ti(i,:)+y1(l,:,i);
            end
        end
        SSw=SSw+ss-1/n1(i)*ti(i,:)'*ti(i,:);
    end
    % Stage I
    %process=[];
    obfun(1:2)=1e150;SSw=trace(SSw);InitialSSw=num2str(SSw)
    cr=.08*SSw/n;
    % Stage II, Contraction-expansion algorithm
    disp(['No.   minSSw'])
    sv=[];
    for v6=1:15
        v5=2;obfun(v5)=obfun(1);rep=0;ne=0;nc=0;
        while v5<=18 && rep<=3
            for v2=0:sqrt(n/2+g*20+v5)
                ti2=zeros(g,k);
                for v=1:12
                    sv(v)=SSw;
                    for p=1:g
                        for q=1:g
                            for l=1:n
                                if q~=p && n1(p)>1
                                    if ~isnan(y1(l,:,p))
                                        ti2(p,:)=ti(p,:)-y1(l,:,p);
                                        ti2(q,:)=ti(q,:)+y1(l,:,p);
                                        dpq=1/(n1(p)-1)*ti2(p,:)*ti2(p,:)'+1/(n1(q)+1)*ti2(q,:)*ti2(q,:)'-1/n1(p)*ti(p,:)*ti(p,:)'-1/n1(q)*ti(q,:)*ti(q,:)';
                                        if dpq>0
                                            ti(p,:)=ti2(p,:);
                                            ti(q,:)=ti2(q,:);
                                            n1(p)=n1(p)-1;
                                            n1(q)=n1(q)+1;
                                            y1(l,:,q)=y1(l,:,p);
                                            y1(l,:,p)=nan(1,k);
                                            SSw=SSw-dpq;
                                            if SSw<obfun(v5)
                                                obfun(v5)=SSw;
                                                ym=y1;
                                                %cagraph
                                                %pause(.00001)
                                            end
                                        end
                                    end
                                end
                            end
                            %cagraph
                            %pause(.00001)
                        end
                    end
                    if sv(v)<=SSw || str2double(str2mat(vpa(obfun(v5),7)))<=str2double(str2mat(vpa(k1ssw(g-1),7)))
                        break
                    end
                    %cagraph
                    %pause(.00001)
                end
                if str2double(str2mat(vpa(obfun(v5),7)))<=str2double(str2mat(vpa(k1ssw(g-1),7)))
                    break
                end
                ti2=zeros(g,k);
                for v=0:2
                    if mod(v,3)==0
                        p1=1;pi=1;p2=g;
                    else
                        p1=g;pi=-1;p2=1;
                    end
                    for p=p1:pi:p2
                        for q=1:g
                            for l=1:n
                                if q~=p && n1(p)>1
                                    if ~isnan(y1(l,:,p))
                                        ti2(p,:)=ti(p,:)-y1(l,:,p);
                                        ti2(q,:)=ti(q,:)+y1(l,:,p);
                                        dpq=1/(n1(p)-1)*ti2(p,:)*ti2(p,:)'+1/(n1(q)+1)*ti2(q,:)*ti2(q,:)'-1/n1(p)*ti(p,:)*ti(p,:)'-1/n1(q)*ti(q,:)*ti(q,:)';
                                        if dpq>-cr
                                            ti(p,:)=ti2(p,:);
                                            ti(q,:)=ti2(q,:);
                                            n1(p)=n1(p)-1;
                                            n1(q)=n1(q)+1;
                                            y1(l,:,q)=y1(l,:,p);
                                            y1(l,:,p)=nan(1,k);
                                            SSw=SSw-dpq;
                                            if SSw<obfun(v5)
                                                obfun(v5)=SSw;
                                                ym=y1;
                                            end
                                            ne=ne+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if mod(v,3)==1
                        cr=cr*exp(-1.7*(ne-.2*n-3*v2-.5*g)/(n+g));
                    end
                end
                if str2double(str2mat(vpa(obfun(v5),7)))<=str2double(str2mat(vpa(k1ssw(g-1),7)))
                    break
                end
                cagraph
                pause(.000001)
                nc=nc+ne;ne=0;
                if mod(v2,2)==0
                    cr=cr*exp(-1.5*(nc-.3*n-5*v2-3*(v5-2)-.75*g)/(n+g));
                else
                    cr=cr*exp(-1.6*(nc-.6*n-8*v2-5*(v5-2)-1.5*g)/(n+g));
                    %disp([nc,cr])
                    nc=0;
                end
                %stage III, Recombination
                %looking for the group with largest SSw
                if rand<.3
                    y1=ym;
                elseif rand>.55
                    ti2=zeros(g,k);
                    for v=1:5
                        sv(v)=SSw;
                        for p=1:g
                            for q=1:g
                                for l=1:n
                                    if q~=p && n1(p)>1
                                        if ~isnan(y1(l,:,p))
                                            ti2(p,:)=ti(p,:)-y1(l,:,p);
                                            ti2(q,:)=ti(q,:)+y1(l,:,p);
                                            dpq=1/(n1(p)-1)*ti2(p,:)*ti2(p,:)'+1/(n1(q)+1)*ti2(q,:)*ti2(q,:)'-1/n1(p)*ti(p,:)*ti(p,:)'-1/n1(q)*ti(q,:)*ti(q,:)';
                                            if dpq>0
                                                ti(p,:)=ti2(p,:);
                                                ti(q,:)=ti2(q,:);
                                                n1(p)=n1(p)-1;
                                                n1(q)=n1(q)+1;
                                                y1(l,:,q)=y1(l,:,p);
                                                y1(l,:,p)=nan(1,k);
                                                SSw=SSw-dpq;
                                                if SSw<obfun(v5)
                                                    obfun(v5)=SSw;
                                                    ym=y1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if sv(v)<=SSw || str2double(str2mat(vpa(obfun(v5),7)))<=str2double(str2mat(vpa(k1ssw(g-1),7)))
                            break
                        end
                    end
                end
                si=zeros(1,g);ti=zeros(g,k);n1=zeros(1,g);
                for i=1:g
                    ss=zeros(k);
                    for l=1:n
                        if ~isnan(y1(l,:,i))
                            ss=ss+x(l,:)'*x(l,:);
                            ti(i,:)=ti(i,:)+x(l,:);
                            n1(i)=n1(i)+1;
                        end
                    end
                    mi(i,:)=ti(i,:)/n1(i);
                    si(i)=trace(ss-1/n1(i)*ti(i,:)'*ti(i,:));
                end
                if g>=7 && mod(v2,2)==1
                    for i=1:g-1
                        for j=i+1:g
                            if si(i)>si(j)
                                tmpD=si(j);
                                si(j)=si(i);
                                si(i)=tmpD;
                                tmpy=y1(:,:,j);
                                y1(:,:,j)=y1(:,:,i);
                                y1(:,:,i)=tmpy;
                            end
                        end
                    end
                    ti=zeros(g,k);n1=zeros(1,g);
                    for i=1:g
                        for l=1:n
                            if ~isnan(y1(l,:,i))
                                ti(i,:)=ti(i,:)+x(l,:);
                                n1(i)=n1(i)+1;
                            end
                        end
                        mi(i,:)=ti(i,:)/n1(i);
                    end
                    midst=1e10;
                    for i=1:g-3
                        for j=i+1:g-2
                            dij(i,j)=(mi(i,:)-mi(j,:))*(mi(i,:)-mi(j,:))';
                            if dij(i,j)<midst
                                midst=dij(i,j);
                                mp1=i;mq1=j;
                            end
                        end
                    end
                    midst=1e10;
                    for i=1:g-3
                        for j=i+1:g-2
                            if dij(i,j)<midst && i~=mp1 && i~=mq1 && j~=mp1 && j~=mq1
                                midst=dij(i,j);
                                mp2=i;mq2=j;
                            end
                        end
                    end
                    % merge the first two groups into one.
                    for l=1:n
                        if ~isnan(y1(l,:,mq1))
                            y1(l,:,mp1)=y1(l,:,mq1);
                            y1(l,:,mq1)=nan(1,k);
                        end
                    end
                    % merge the second two groups into one.
                    for l=1:n
                        if ~isnan(y1(l,:,mq2))
                            y1(l,:,mp2)=y1(l,:,mq2);
                            y1(l,:,mq2)=nan(1,k);
                        end
                    end
                    %divide the first largest group into two groups
                    l=0;
                    while l<=n
                        l=random('unid',n,1);
                        if ~isnan(y1(l,:,g))
                            y1(l,:,mq1)=y1(l,:,g);
                            y1(l,:,g)=nan(1,k);
                            l=n+1;
                        end
                    end
                    %divide the second largest group into two groups
                    l=0;
                    while l<=n
                        l=random('unid',n,1);
                        if ~isnan(y1(l,:,g-1))
                            y1(l,:,mq2)=y1(l,:,g-1);
                            y1(l,:,g-1)=nan(1,k);
                            l=n+1;
                        end
                    end
                elseif g>=3
                    if rand<.7 % searching for the largest group (max SSw(i)).
                        mxd=0;
                        for i=1:g
                            if si(i)>mxd
                                mxd=si(i);
                                mxq=i;
                            end
                        end
                    else % searching for the largest group (max n1(i)).
                        mxd=0;
                        for i=1:g
                            if n1(i)>mxd
                                mxd=n1(i);
                                mxq=i;
                            end
                        end
                    end
                    if rand<=.7
                        % looking for the groups with minimal distance.
                        midst=1e10;
                        for i=1:g-1
                            for j=i+1:g
                                di=(mi(i,:)-mi(j,:))*(mi(i,:)-mi(j,:))';
                                if di<midst && i~=mxq && j~=mxq
                                    midst=di;
                                    mp=i;mq=j;
                                end
                            end
                        end
                    else % the two merging group is determined randomly.
                        mp=random('unid',g,1);
                        while mp==mxq
                            mp=random('unid',g,1);
                        end
                        mq=random('unid',g,1);
                        while mq==mp || mq==mxq
                            mq=random('unid',g,1);
                        end
                    end
                    % merge two groups into one.
                    for l=1:n
                        if ~isnan(y1(l,:,mq))
                            y1(l,:,mp)=y1(l,:,mq);
                            y1(l,:,mq)=nan(1,k);
                        end
                    end
                    %divide the largest group into two groups
                    l=0;
                    while l<=n
                        l=random('unid',n,1);
                        if ~isnan(y1(l,:,mxq))
                            y1(l,:,mq)=y1(l,:,mxq);
                            y1(l,:,mxq)=nan(1,k);
                            l=n+1;
                        end
                    end
                end
                %recalculating the SSw
                SSw=zeros(k);ti=zeros(g,k);n1=zeros(1,g);
                for i=1:g
                    ss=zeros(k);
                    for l=1:n
                        if ~isnan(y1(l,:,i))
                            ss=ss+x(l,:)'*x(l,:);
                            ti(i,:)=ti(i,:)+x(l,:);
                            n1(i)=n1(i)+1;
                        end
                    end
                    SSw=SSw+ss-1/n1(i)*ti(i,:)'*ti(i,:);
                end
                SSw=trace(SSw);
            end
            if str2double(str2mat(vpa(obfun(v5),7)))==str2double(str2mat(vpa(obfun(v5-1),7)));
                rep=rep+1;
            else
                rep=0;
            end
            cagraph
            pause(.000001)
            v5=v5+1;obfun(v5)=obfun(v5-1);
            cr=(.15*obfun(v5)/(n+g)+cr/5)/2;ne=0;nc=0;
            disp([num2str(v5-2),'    ',num2str(obfun(v5))])
            if str2double(str2mat(vpa(obfun(v5),7)))<=str2double(str2mat(vpa(k1ssw(g-1),7)))
                break
            end
        end % next v5
        wt=input('Another round of standardized clustering? (Yes)=1/(No)=0 ');
        %wt=0;
        if wt==0
            break
        end
        ti=zeros(g,k);n1=zeros(1,g);SSw=zeros(k);
        for i=1:g
            ss=zeros(k);
            for l=1:n
                if ~isnan(ym(l,:,i))
                    ss=ss+x1(l,:)'*x1(l,:);
                    ti(i,:)=ti(i,:)+x1(l,:);
                    n1(i)=n1(i)+1;
                end
            end
            SSw=SSw+ss-1/n1(i)*ti(i,:)'*ti(i,:);
        end
        SSw=SSw/(n-g);
        jj=ones(n,1);sd=sqrt(diag(SSw)');
        for j=1:k
            if sd(j)<=0
                sd(j)=0.001;
            end
        end
        jj=jj*sd;
        x=x1./jj;
        SSw=zeros(k);y1=nan(n,k,g);ti=zeros(g,k);n1=zeros(1,g);
        for i=1:g
            ss=zeros(k);
            for l=1:n
                if ~isnan(ym(l,:,i))
                    y1(l,:,i)=x(l,:);
                    ss=ss+x(l,:)'*x(l,:);
                    ti(i,:)=ti(i,:)+x(l,:);
                    n1(i)=n1(i)+1;
                end
            end
            SSw=SSw+ss-1/n1(i)*ti(i,:)'*ti(i,:);
        end
        SSw=trace(SSw);cr=.003*SSw/(n+g);ym=y1;obfun(1:2)=SSw+1e10;
    end % next v6
    SSw=0;ti=zeros(g,k);n1=zeros(1,g);
    for i=1:g
        ss=0;s(i)=0;
        for l=1:n
            if ~isnan(ym(l,:,i))
                ss=ss+x1(l,:)'*x1(l,:);
                ti(i,:)=ti(i,:)+x1(l,:);
                n1(i)=n1(i)+1;
            end
        end
        s(i)=trace(ss-1/n1(i)*ti(i,:)'*ti(i,:));
        SSw=SSw+s(i);
    end
    cagraph
    pause(.001)
    minSSw=trace(SSw);
    TDw=0;TSDw=0;
    for i=1:g
        for l1=1:n-1
            for l2=l1+1:n
                if ~isnan(ym(l1,:,i)) & ~isnan(ym(l2,:,i))
                    TDw=TDw+sqrt((ym(l1,:,i)-ym(l2,:,i))*(ym(l1,:,i)-ym(l2,:,i))');
                    TSDw=TSDw+(ym(l1,:,i)-ym(l2,:,i))*(ym(l1,:,i)-ym(l2,:,i))';
                end
            end
        end
    end
    indivals=[];% individuals in each group
    for i=1:g
        m=1;
        for l=1:n
            if sum(abs(ym(l,:,i)))>0
                m=m+1;
                indivals(i,m)=l;
            end
        end
    end
    disp('Clustering results:')
    disp('group     individuals...')
    for i=1:g
        ind=indivals(i,:);ind(ind==0)=[];
        disp([num2str(i),': ',num2str(ind)])
    end
    for i=1:g
        mg(i,:)=ti(i,:)./n1(i);
    end
    groupmean=[];
    for i=1:g
        groupmean=[groupmean; i,n1(i),mg(i,:),s(i)];
    end
    len=size(groupmean,2);
    groupmean=vpa(groupmean,7);
    traits=[' Group ';'Number '];
    for j=1:k
        if j<10
            traits=[traits;' X',num2str(j),'mean'];
        else
            traits=[traits;'X',num2str(j),'mean'];
        end
    end
    traits=[traits;' SSw(i)'];
    line1=vpa(ones(1,len),6);
    groupmean=[line1;groupmean];
    groupmean(1,1)=traits(1,:);
    groupmean(1,2)=traits(2,:);
    for i=1:k
        groupmean(1,i+2)=traits(i+2,:);
    end
    groupmean(1,k+3)=str2sym(traits(k+3,:));
    groupmean(g+2,1)='Total ';groupmean(g+2,2)=vpa(n,6);
    tm=mean(x1); % total mean
    for i=1:k
        groupmean(g+2,i+2)=vpa(tm(i),6);
    end
    groupmean(g+2,k+3)=vpa(minSSw,6)
    Final_SSw____TDw____TSDw=[minSSw,TDw,TSDw]
    OF(g,init)=obfun(v5);tm(g,init)=cputime-ta1;
end