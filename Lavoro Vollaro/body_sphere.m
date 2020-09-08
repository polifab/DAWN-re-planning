function body_sphere(obj_id,obj_pos)
% BODY_SPHERE(obj_id,rr) plots a sphere at coordinates OBJ_POS
%   and applies the image of OBJ_ID body to its surface.
%   Available: all solar planets, Pluto, Vesta, Ceres, Sun.
%
%   obj_id   - numeric identifier of the target body (1-12)
%   obj_pos  - (x,y,z) coordinates of target body, wrt the Sun
%

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
        
    radii = [2439.5
             6052 
             6378
             3396
             71492
             60268
             25559
             24764
             1185
             262.7
             476.2
             695508]; %[km]
         
	R = radii(obj_id); %[km]
    
    [xx,yy,zz] = sphere(100);
    sp_hand = surface(obj_pos(1)+R*xx,obj_pos(2)+R*yy,obj_pos(3)+R*zz);
    
    img = imread(body(obj_id));
    set(sp_hand,'facecolor','texture',...
        'cdata',im2double(imrotate(img,180)),...
        'edgecolor','none');
end