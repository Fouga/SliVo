%saves txt file with feature's data to output_file


function [des1,loc1]=save_keys(output_file,Im_1, method,thresh)


% p = inputParser;
% p.addRequired('output_file');
% p.addRequired('image');
% 
% defaultMethod = 'SURF';
% defaultThresh = 0.0002; 
% 
% p.addParameter('thresh', defaultThresh, @(x) isnumeric(x) && length(x)==1);
% p.addOptional('meth', defaultMethod, @(x) ischar(x));
% 
% parse(p, output_file,image, varargin{:}); %Pass the values of all of the function inputs.
% thresh = p.Results.thresh;% access the valuse of threshold
% % method = p.Results.meth;
% if nargin < 4
%   thresh = 0.0002;
% else thresh = varargin{1}
% end
  

if (strcmp(method, 'SIFT'))
    cd ('./sift_dent_project');
    [im1, des1, loc1] = sift2(Im_1);
    copyfile('tmp.key', output_file); %copy file with keys to the new dir
    
elseif (strcmp(method,'ASIFT'))
    cd ('./demo_ASIFT_src');
    [im1]=asift_out(Im_1);
    copyfile('keys1.txt', output_file);
    
elseif (strcmp(method,'SURF'))
    
    Options.tresh=thresh; %Hessian response treshold (default 0.0002)
    Ipts1=OpenSurf(Im_1,Options); %key points of the image
    
    des1_surf = reshape([Ipts1.descriptor],64,[]); 
    des1=des1_surf'; % matrix K by 64
    loc1.y=reshape([Ipts1.y],1,[])';% x = col, y - row
    loc1.x=reshape([Ipts1.x],1,[])'; %locs of  [column, row]   
    scale=reshape([Ipts1.scale],1,[])';
    orientation=reshape([Ipts1.orientation],1,[])';
    keys_number=size(des1,1);
    keys_data=zeros(keys_number,size(des1,2)+4); %matrix cols=64+4, rows=amount of detected features
    keys_data(1:keys_number,1:4)=[loc1.x, loc1.y,scale, orientation];
    keys_data(1:keys_number,5:end)=des1;

    header=[keys_number, size(des1,2)];
 
    dlmwrite(output_file,header,'delimiter',' ');
    dlmwrite(output_file,keys_data,'-append','delimiter',' ');
   
end

end




