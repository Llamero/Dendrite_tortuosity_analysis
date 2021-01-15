tortuosity = 1; //(0 = straight line, 1 = random line)
n_lines = 240; //number of lines to draw
segment_length = 5; //Length of each line segment
min_segments = 5; //Number of line segments to stop drawing line if not yet off edge of image
max_segments = 100;
image_size = 1024;
random("seed", 10); //Random seed
if(!isOpen("Test image")) newImage("Test image", "8-bit black", image_size, image_size, 2);
setColor(255);

for(a=0; a<n_lines; a++){
	x1 = round(random()*image_size)-1;
	y1 = round(random()*image_size)-1;
	n_segments = round(random()*(max_segments-min_segments))+min_segments;
	for(b=0; b<n_segments; b++){
		bias_radians = atan2(y1-(image_size/2), x1-(image_size/2))-PI/2;
		rad = bias_radians + ((random()-0.5)*2*PI*tortuosity);
		x2 = x1 + round(segment_length*cos(rad));
		y2 = y1 + round(segment_length*sin(rad));
		drawLine(x1, y1, x2, y2);
		x1 = x2;
		y1 = y2;

		//Wrap around frame if line hits edge
		if(x1 < 0) x1 = image_size -1;
		if(y1 < 0) y1 = image_size -1;
		if(x1 >= image_size) x1 = 0;
		if(y1 >= image_size) y1 = 0;
	}
}
title = "Dendrite mask - Spiral " + tortuosity;
title = replace(title, "\\.", "_");
rename(title);
