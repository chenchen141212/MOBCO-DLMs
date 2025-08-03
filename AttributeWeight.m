
function w=AttributeWeight(data,kmedoidsclusters,kmedoidscenter_)

[~,m]=size(data);
Dj=zeros(1,m);
d=zeros(1,m);
a=8;
    for j=1:m
        for k=1:size(kmedoidscenter_,1)        
            xj=data(kmedoidsclusters==k,j);
            yj=kmedoidscenter_(k,j);
            djk=sum(xj==yj);
        end
        Dj(1,j)=sum(djk);
    end
%    w=Dj./sum(Dj);

     for j=1:m
         if Dj(1,j)~=0 
            for z=1:m 
                if Dj(1,z)~=0
                    delta=1/(a-1);
                    d(1,z)=(Dj(1,j)/Dj(1,z))^delta; 
                else
%                    do noting
                end              
            end
            w(1,j)=1/sum(d);
         else
            w(1,j)=0;
         end
     end
end 