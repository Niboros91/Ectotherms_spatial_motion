function [Z,Z_wind,Link]=Adjency_matrix(Columns,Rows)
%We use as input the number points per column and row
N=Columns*Rows;
Pos=zeros(N,2);
pos=1;

%To know file and column
for i=1:Rows
    for j=1:Columns
        Pos(pos,1) = j;
        Pos(pos,2) =i;
        pos=pos+1; 
    end
end
Z=zeros(N,N); %We don't add any extra point
Z_wind = zeros(N,N);
Link=zeros(N,2);
e=1;  

%To stablish which squares ae neighbors and which not
for k=1:N
    for i=1:N
     
            %To  check if they are neighbours
            if (Pos(k,1)==Pos(i,1)) 
                if (abs(Pos(k,2)-Pos(i,2))<2) && (Pos(k,2)~=Pos(i,2))
                    Z(k,i)=1;
                    Link(e,1)=k;
                    Link(e,2)=i;
                    e=e+1;
                    % x coordinate equal
                    if (Pos(k,2)-Pos(i,2)) > 0
                        Z_wind(k,i) = - 1;
                    else
                        Z_wind(k,i) = 1;
                    end
                    
                end
            end
             if (Pos(k,2)==Pos(i,2)) 
                 if (abs(Pos(k,1)-Pos(i,1))<2) && (Pos(k,1)~=Pos(i,1))
                   Z(k,i)=1; 
                   Link(e,1)=k;
                   Link(e,2)=i;
                   e=e+1;
                   
                   % y coordinate equal
                    if (Pos(k,1)-Pos(i,1)) > 0
                        Z_wind(k,i) = - 2;
                    else
                        Z_wind(k,i) = 2;
                    end
                   
                 end
             end
        
    end
end


end
