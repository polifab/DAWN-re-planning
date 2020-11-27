function body_sphere(obj_id,obj_pos)
% BODY_SPHERE(obj_id,rr) plots a sphere at coordinates OBJ_POS
%   and applies the image of OBJ_ID body to its surface.
%   Available: all solar planets, Pluto, Vesta, Ceres, Sun.
%
%   obj_id   - numeric identifier of the target body (1-12)
%
%   obj_pos  - (x,y,z) coordinates of target body, wrt the Sun
%

    %% Constants
    body = ["mercury.jpg" 
            "venus.jpg" 
            "earth.jpg" 
            "mars.jpg" 
            "jupiter.jpg" 
            "saturn.jpg" 
            "uranus.jpg" 
            "neptune.jpg" 
            "pluto.jpg" 
            "vesta.jpg" 
            "ceres.jpg" 
            "sun.jpg"];
        
    radii = [2439.7
             6051.8 
             6371
             3389.5
             69911
             58232
             25362
             24622
             1151
             262.7
             476.2
             695508]; %[km]
         
	R = radii(obj_id); %[km]
    
    %% Sphere creation
    [xx,yy,zz] = sphere(100);
    sp_hand = surface(obj_pos(1)+R*xx,obj_pos(2)+R*yy,obj_pos(3)+R*zz);
    
    %% Surface change
    img = imread(body(obj_id));

	%img = imrotate(img,180);
% 	tform = affine2d([1 0 0; 0 -1 0; 0 0 1]);
	tform = affine3d(...
			[[1 0 0; 0 -1 0; 0 0 1], [0 0 0]'; 0 0 0 1] );
		
	img = imwarp(img,tform);
	
    set(sp_hand,'facecolor','texture',...
        'cdata',im2double(img),...
        'edgecolor','none');
    
end