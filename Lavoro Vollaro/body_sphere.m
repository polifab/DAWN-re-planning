function body_sphere(obj_id,rr)

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
         
	R = radii(obj_id);
    
    [xx,yy,zz] = sphere(100);
    sp_hand = surface(rr(1)+R*xx,rr(2)+R*yy,rr(3)+R*zz);
    
    img = imread(body(obj_id));
    % Set it on SUN
    set(sp_hand,'facecolor','texture',...
        'cdata',im2double(imrotate(img,180)),...
        'edgecolor','none');
end