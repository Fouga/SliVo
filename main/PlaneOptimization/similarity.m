function out = similarity(A, B, method)

[m n]=size(A);

A = round(A); A(find(A<0))=0; A(find(A>255))=255;
B = round(B); B(find(B<0))=0; B(find(B>255))=255;
switch method
   case 'MSD'
        out = sum(sum((A-B).^2))/(m*n);

    case 'MAD'
		% Task 1
		out = sum(sum(abs(A-B)))/(m*n);

    case 'CORR'
		% Task 2
		out =  sum(sum(A.*B))/(sum(sum(A)));

    case 'CC'
		% Task 3
        muA=mean(mean(A));
        muB=mean(mean(B));
		out = sum(sum((A-muA).*(B-muB)))/sqrt(sum(sum((A-muA).^2.*(B-muB).^2)));
        
    case 'JE'
		% Task 4
		a=A(:)+1; %grey values
        b=B(:)+1;
        Jhist=zeros(256,256);
        for i=1:m*n
            Jhist(a(i),b(i))=Jhist(a(i),b(i))+1;
        end
        Pr = Jhist./(m*n);
        out = -sum(sum(Pr.*(log2(Pr+(Pr==0)))));  % joint entropy
        % (Pr+(Pr==0)) returns 1 if Pr=0, so log is determined 
    
    case 'NMI'
% Task 5
        a=A(:)+1; %grey values
        b=B(:)+1;
        Jhist=zeros(256,256);
        for i=1:m*n
            Jhist(a(i),b(i))=Jhist(a(i),b(i))+1;
        end
        %figure, imshow(log(Jhist),[]); % x - A, y - B
        Pr = Jhist./(m*n);
        %entropy of the first image
        EntA=0;PrA=sum(Pr,1);
        for i=1:256    %  column 
            if( PrA(i)==0 )
            else
            EntA = EntA -(PrA(i)*(log2(PrA(i)))); %marginal entropy for image 1
            end
        end
        EntB=0;PrB=sum(Pr,2);
        for j=1:256;    %rows
            if( PrB(j)==0 )
            else
                EntB = EntB + -(PrB(j)*(log2(PrB(j)))); %marginal entropy for image 2
            end   
        end
        
        Ent = -sum(sum(Pr.*(log2(Pr+(Pr==0))))); % joint entropy
       
        
        out = (EntA + EntB)/ Ent; % Mutual information
        
 
    case 'MIND'
        
        mind2 = MIND_descriptor2D(A,1);
        mind1 = MIND_descriptor2D(B,1);
        out = -sum(sum(sum(abs(mind1 - mind2),3)))/1000;
    case 'NCC'
        
        C =  normxcorr2_general(A,B);
        out = max(C(:));
        
end