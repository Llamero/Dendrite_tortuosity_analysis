setBatchMode(true);
close("*");
n_steps = 100;
n_lines = 240; //number of lines to draw
segment_length = 5; //Length of each line segment
min_segments = 20; //Number of line segments to stop drawing line if not yet off edge of image
max_segments = 100;
image_width = 1024;
image_height = 1024;
bias_angle = 23; //Predominant direction of neurites
bias_radians = bias_angle*PI/180;

//Initialize global variables
tolerance = -1;
size = -1;
label_array = newArray("Mean", "Median");

//Direction variables
//FFT_substack_nSlices = 3; //Thickness in slices of substack to generate FFT (ensures neurites on different slices are included);
//n_stdev = 2; //Number of standard deviations above the mean to threshold the FFT
tolerance_range = 0.1; //+/- factor range to test tolerance
tolerance_step = tolerance_range; //Step size to test tolerance
min_size_cutoff = 4; //Smallest size to be exluded in particle analyzer
max_size_cutoff = 20; //Largest size to be exluded in particle analyzer
FFT_min_radius = 1/64; //Ratio of total width of FFT to set as min cutoff - i.e. lowest sptatial frequency
FFT_max_radius = 1/4; //Ratio of total width of FFT to set as max cutoff - i.e. highest sptatial frequency
sector_width = 10; //Width of sector off FFT to look for directionality
n_degrees = 360; //Number of degrees to sample in FFT
grid_division = 4; //Needs to be a power of 2

i = simulateSpiral();
//i = simulateLine();
makeRectangle(image_width*0.375, image_width*0.75, image_width*0.25, image_width*0.25);
run("Crop");
//createFFT(i);
getDirectionality(i);

function simulateLine(){
	newImage("Test image", "8-bit black", image_width, image_height, n_steps+1);
	setColor(255);
	slice = 1;
	for(tortuosity=0; tortuosity<=1; tortuosity += 1/n_steps){
		showProgress(slice, n_steps); 
		random("seed", 10); //Random seed
		Stack.setSlice(slice++);
		setMetadata("Label", "Toruosity = " + tortuosity);
		for(a=0; a<n_lines; a++){ //Draw random starting point and number of line segments
			x1 = round(random()*image_width)-1;
			y1 = round(random()*image_height)-1;
			n_segments = round(random()*(max_segments-min_segments))+min_segments;
			for(b=0; b<n_segments; b++){ //Draw random angle for each segment
				rad = bias_radians + ((random()-0.5)*2*PI*tortuosity);
				x2 = x1 + round(segment_length*cos(rad));
				y2 = y1 + round(segment_length*sin(rad));
				drawLine(x1, y1, x2, y2);
				x1 = x2;
				y1 = y2;
		
				//Wrap around frame if line hits edge
				if(x1 < 0) x1 = image_width -1;
				if(y1 < 0) y1 = image_height -1;
				if(x1 >= image_width) x1 = 0;
				if(y1 >= image_height) y1 = 0;
			}
		}
	}
	title = "Dendrite mask - Linear " + tortuosity;
	title = replace(title, "\\.", "_");
	rename(title);
	return title;
}

function simulateSpiral(){
	newImage("Test image", "8-bit black", image_width, image_height, n_steps);
	setColor(255);
	slice = 1;
	for(tortuosity=0; tortuosity<=1; tortuosity += 1/n_steps){
		showProgress(slice, n_steps); 
		random("seed", 10); //Random seed
		Stack.setSlice(slice++);
		setMetadata("Label", "Toruosity = " + tortuosity);
		for(a=0; a<n_lines; a++){
			x1 = round(random()*image_width)-1;
			y1 = round(random()*image_height)-1;
			n_segments = round(random()*(max_segments-min_segments))+min_segments;
			for(b=0; b<n_segments; b++){
				bias_radians = atan2(y1-(image_height/2), x1-(image_width/2))-PI/2;
				rad = bias_radians + ((random()-0.5)*2*PI*tortuosity);
				x2 = x1 + round(segment_length*cos(rad));
				y2 = y1 + round(segment_length*sin(rad));
				drawLine(x1, y1, x2, y2);
				x1 = x2;
				y1 = y2;
		
				//Wrap around frame if line hits edge
				if(x1 < 0) x1 = image_width -1;
				if(y1 < 0) y1 = image_height -1;
				if(x1 >= image_width) x1 = 0;
				if(y1 >= image_height) y1 = 0;
			}
		}
	}
	title = "Dendrite mask - Spiral " + tortuosity;
	title = replace(title, "\\.", "_");
	rename(title);
	return title;
}



function createFFT(i){
	run("FFT Options...", "raw");
	selectWindow(i);
	slices = nSlices;
	for(a=1; a<=slices; a++){
		selectWindow(i);
		setSlice(a);
		run("FFT");
		selectWindow("PS of " + i);
		run("Log");
		if(a==1){
			rename("FFT Stack");
		}
		else{
			run("Concatenate...", "  title=[FFT Stack] image1=[FFT Stack] image2=[PS of " + i + "] image3=[-- None --]");
		}
	}
}

function getDirectionality(image){
	run("FFT Options...", "raw");
	selectWindow(i);
	slices = nSlices;
	newImage("Statistics", "32-bit black", n_degrees, (2*label_array.length)+1, 1);
	for(slice=1; slice<=slices-1; slice++){
		print(slice + " of " + slices);
		selectWindow(i);
		setSlice(slice);
		run("FFT");
		selectWindow("PS of " + i);
		rename("FFT1");
		//run("Log");
		//run("Gaussian Blur...", "sigma=8");
	
		//make angle map
		selectWindow("FFT1");
		if(slice == 1){
			getDimensions(FFT_width, dummy, dummy, dummy, dummy);
			newImage("Angle", "8-bit black", FFT_width, FFT_width, 1);
			selectWindow("Angle");
			setPixel(FFT_width/2, FFT_width/2, 255);
			run("Invert");
			run("Exact Euclidean Distance Transform (3D)");
			selectWindow("EDT");
			run("Canny Edge", "xsize=1 ysize=1 zsize=1 x-sigma,=1 y-sigma,=1 z-sigma,=1 low=0.2000000 high=0.6000000 angle scaling=-1");
		
			//Create mask for fixed annular ROI to find direction preference
			selectWindow("EDT");
			setThreshold(FFT_width*FFT_min_radius, FFT_width*FFT_max_radius);
			run("NaN Background");
			imageCalculator("Divide", "EDT","EDT");
	
			//Set angle map to be scaled in degrees and apply distance ROI
			selectWindow("canny phi");
			run("Rotate 90 Degrees Left");
			run("Multiply...", "value=1.41176470588");
			imageCalculator("Multiply", "canny phi","EDT");
			close("EDT");
			close("Angle");
		}
	
		//Measure statistics for each sector of the FFT
		selectWindow("FFT1");
		run("Select None");
		for(b=0; b<4; b++){
			for(a=45; a<135; a++){
				angle = (a+b*90)%360;
				selectWindow("canny phi");
				run("Select None");
				//Measure for 45 to 135 - then tranform angle map 90° and repeat 4x
				setThreshold(a-sector_width/2, a+sector_width/2);
				run("Create Selection");
				selectWindow("FFT1");
				run("Restore Selection");
				List.setMeasurements();
				run("Select None");
				selectWindow("Statistics");
				for(c=0; c<label_array.length; c++){
					setPixel(angle, c, parseFloat(List.get(label_array[c])));
				}
			}
			selectWindow("canny phi");	
			run("Rotate 90 Degrees Left");
		}
	
		//Output results to results table
		newImage("Polar plots", "RGB black", 1024, 1024, label_array.length);
		setForegroundColor(255, 0, 255);
		makeRectangle(511, 0, 1, 1024);
		run("Fill", "stack");
		makeRectangle(0, 511, 1024, 1);
		run("Fill", "stack");
		run("Select None");
		selectWindow("Statistics");
		results_array = newArray(n_degrees);
		line_array = newArray(4);
		for(a=1; a<label_array.length; a++){
			//Measure FFT aspect ratios
			selectWindow("Statistics");
			makeRectangle(0, a, n_degrees, 1);
			getStatistics(dummy, dummy, min, max);
	
			//Create polar plot of FFT aspect ratios
			scale_radius = 512/max;
			selectWindow("Statistics");
			line_array[0] = getPixel((n_degrees-1),a);
			line_array[1] = getPixel((n_degrees-1),a);
			line_array[0] = cos((n_degrees-1-90)*PI/180)*line_array[0]*scale_radius+511;
			line_array[1] = 511-sin((n_degrees-1-90)*PI/180)*line_array[1]*scale_radius; //positive Y is pointing down
			setForegroundColor(0, 255, 0);
			selectWindow("Polar plots");
			Stack.setSlice(a+1);
			setMetadata("Label", label_array[a]);
			max_ratio = -1;
			max_angle = -1;
			for(b=0; b<n_degrees; b++){
				selectWindow("Statistics");
				ratio = getPixel(b, a)/getPixel((b+90)%n_degrees, a); //Measure ratio of orthogonal vectors
				if(ratio > max_ratio){
					max_ratio = ratio;
					max_angle = (b+90)%n_degrees;
				}
				line_array[2] = getPixel(b,a);
				line_array[3] = getPixel(b,a);
				line_array[2] = cos((b-90)*PI/180)*line_array[2]*scale_radius+511;
				line_array[3] = 511-sin((b-90)*PI/180)*line_array[3]*scale_radius;
				selectWindow("Polar plots");
				makeLine(line_array[0], line_array[1], line_array[2], line_array[3]);
				run("Draw", "slice"); 
				line_array[0] = line_array[2];
				line_array[1] = line_array[3]; 			
			}
	
			//Save results to output matrix
			selectWindow("Statistics");
			setPixel(0, label_array.length + a, max_ratio);
			setPixel(1, label_array.length + a, max_angle);
	
			//Get polar plot statistics
			selectWindow("Polar plots");
			run("Duplicate...", "title=plot");
			selectWindow("plot");
			run("Split Channels");
			close("plot (red)");
			close("plot (blue)");
			selectWindow("plot (green)");
			run("Invert");
			run("Dilate");
			run("Fill Holes");
			run("Erode");
			run("Create Selection");
			List.setMeasurements;
			ellipse_ratio = parseFloat(List.get("Major"))/parseFloat(List.get("Minor"));
			ellipse_angle = PI/180*parseFloat(List.get("Angle"));
			ellipse_major = parseFloat(List.get("Major"));
			ellipse_minor = parseFloat(List.get("Minor"));
			close("plot (green)");
			selectWindow("Polar plots");
			dx = cos(ellipse_angle)*ellipse_major/2;
			dy = sin(-1*ellipse_angle)*ellipse_major/2;
			
			//Draw fitted ellipse
			newImage("Ellipse plots", "RGB black", 1024, 1024, 1);
			setMetadata("Label", "Tortuosity = " + 1/ellipse_ratio);
			setForegroundColor(255, 0, 255);
			makeRectangle(511, 0, 1, 1024);
			run("Fill", "stack");
			makeRectangle(0, 511, 1024, 1);
			run("Fill", "stack");
			run("Select None");
			setForegroundColor(0, 255, 0);
			makeEllipse(511-dx, 511-dy, 511+dx, 511+dy, 1/ellipse_ratio);
			run("Draw", "slice");
			run("Select None");

			selectWindow("Polar plots");
			setSlice(1);
			run("Delete Slice");
			if(slice == 1){
				selectWindow("Polar plots");
				rename("Polar stack");
				selectWindow("Ellipse plots");
				rename("Ellipse stack");	
				selectWindow("FFT1");
				rename("FFT stack");	
			}
			else{
				run("Concatenate...", "  title=[Polar stack] image1=[Polar stack] image2=[Polar plots] image3=[-- None --]");
				run("Concatenate...", "  title=[Ellipse stack] image1=[Ellipse stack] image2=[Ellipse plots] image3=[-- None --]");
				run("Concatenate...", "  title=[FFT stack] image1=[FFT stack] image2=[FFT1] image3=[-- None --]");
			}
		}
//		close("Polar plots");
//		close("FFT1");
		selectWindow("Statistics");
		run("Select None");
		selectWindow("canny phi");
		run("Select None");
	}
	close("canny phi");
	close("Statistics");
}
setBatchMode("exit and display");
