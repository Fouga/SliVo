matlabpool('open',4);
siz=181;
sigma = siz/2.34/2;
siz = (siz-1)./2;
x = linspace(-siz, siz, siz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
matr_smoothX = convn(matr,gaussFilter,'same');
matr_smoothY = convn(matr_smoothX,gaussFilter','same');

siz=3;
%sigma = 0.65;
sigma = siz/2.34/2;
siz = (siz-1)./2;
x = linspace(-siz, siz, siz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); 
z=zeros(1,1,numel(x));
z(:,:,:)=gaussFilter; 
matr_smoothZ = convn(matr_smoothY,z,'same');

clear matr_smoothX
clear matr_smoothY


figure,imshow(matr_smoothZ(:,:,209),[])
