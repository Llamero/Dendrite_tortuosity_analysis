close("*");

//Initialize global variables
tolerance = -1;
size = -1;
label_array = newArray("Mean", "Median", "Min", "Max");

//Direction variables
//FFT_substack_nSlices = 3; //Thickness in slices of substack to generate FFT (ensures neurites on different slices are included);
//n_stdev = 2; //Number of standard deviations above the mean to threshold the FFT
tolerance_range = 0.1; //+/- factor range to test tolerance
tolerance_step = tolerance_range; //Step size to test tolerance
min_size_cutoff = 4; //Smallest size to be exluded in particle analyzer
max_size_cutoff = 20; //Largest size to be exluded in particle analyzer
FFT_min_radius = 1/128; //Ratio of total width of FFT to set as min cutoff - i.e. lowest sptatial frequency
FFT_max_radius = 1/32; //Ratio of total width of FFT to set as max cutoff - i.e. highest sptatial frequency
sector_width = 360; //Width of sector off FFT to look for directionality
n_degrees = 360; //Number of degrees to sample in FFT
grid_division = 1; //Needs to be a power of 2

setBatchMode(true);
in_dir = getDirectory("INPUT directory:");
file_list = getFileList(in_dir);
out_dir = getDirectory("OUTPUT directory:");

for(a=0; a<file_list.length; a++){
	if(startsWith(file_list[a], "Dendrite mask - ")){
		open(in_dir + file_list[a]);
		image = getTitle();
		print("Processing image: " + image + ".  " + a+1 + " of " + file_list.length+1 + " files.");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("MAX_" + image);
		rename("MAX_Dendrite mask");
		getDimensions(width, dummy, dummy, dummy, dummy);
		rect_size = width/grid_division;
		rect_step = rect_size;
		newImage("Quiver Vectors", "32-bit grayscale-mode label", grid_division, grid_division, 3, label_array.length, 1);
		newImage("Statistics", "32-bit black", n_degrees, (2*label_array.length)+1, 1);
		total_progress = grid_division*grid_division;
		for(x=0; x<grid_division; x++){
			for(y=0; y<grid_division; y++){
				current_progress = x*grid_division + y + 1;
				print("Processing grid " + current_progress + " of " + total_progress);
				selectWindow("MAX_Dendrite mask");
				run("Grays"); //Nead for correct measurement of mean
				makeRectangle(x*rect_step, y*rect_step, rect_size, rect_size);

				//Measure neurite area
				getStatistics(dummy, neurite_mean);
				neurite_area = neurite_mean/255;
				neurite_area *= rect_size*rect_size;
				getDirectionality(image, x, y, neurite_area);
			}
		}
		selectWindow("Quiver Vectors");
		image = replace(image, "\\..+?$", ""); 
		image = replace(image, "Dendrite mask - ", ""); 
		saveAs("tiff", out_dir + "Quiver Vectors - " + image);
		run("Duplicate...", "title=Ratio duplicate channels=1 slices=2");
		selectWindow("Quiver Vectors - " + image + ".tif");
		run("Duplicate...", "title=Area duplicate channels=3 slices=2");
		selectWindow("Ratio");
		run("Select All");
		List.setMeasurements;
		mean = List.getValue("Mean");
		median = List.getValue("Median");
		setResult("Sample", nResults, image);
		setResult("Mean", nResults-1, mean);
		setResult("Median", nResults-1, median);
		selectWindow("Area");
		getStatistics(image_size, mean_area);
		total_area = image_size * mean_area;
		run("Divide...", "value=" + total_area); //Normalize area to 1
		imageCalculator("Multiply", "Ratio","Area");
		selectWindow("Ratio");
		run("Bin...", "x=" + grid_division + " y=" + grid_division + " bin=Sum");
		weighted_ratio = getPixel(0,0);
		setResult("Weighted", nResults-1, weighted_ratio);
		updateResults();
		close("*");
	}
}

function getDirectionality(image, x_pos, y_pos, neurite_area){
	//Get FFT of mask
	run("FFT Options...", "raw do");
	selectWindow("MAX_Dendrite mask");
	run("FFT");
	selectWindow("PS of MAX_Dendrite mask");
	rename("FFT1");
	//run("Log");
	//run("Gaussian Blur...", "sigma=8");

	//make angle map
	selectWindow("FFT1");
	if(x_pos == 0 && y_pos == 0){
		getDimensions(FFT_width, dummy, dummy, dummy, dummy);
		newImage("Angle", "8-bit black", FFT_width, FFT_width, 1);
		selectWindow("Angle");
		setPixel(FFT_width/2, FFT_width/2, 255);
		run("Invert");
		run("Exact Euclidean Distance Transform (3D)");
		selectWindow("EDT");
		run("Canny Edge", "xsize=1 ysize=1 zsize=1 x-sigma,=1 y-sigma,=1 z-sigma,=1 low=0.2000000 high=0.6000000 angle scaling=-1");
	
		//Create mask for fixed annualr ROI to find direction preference
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
	for(a=0; a<label_array.length; a++){
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
		run("Fill Holes");
		run("Create Selection");

		List.setMeasurements;
		ellipse_ratio = parseFloat(List.get("Major"))/parseFloat(List.get("Minor")); 
		selectWindow("Quiver Vectors");
		Stack.setSlice(a+1);
		Stack.setChannel(1);
		setMetadata("Label", label_array[a] + " - Ratio");
		setPixel(x, y, ellipse_ratio);
		Stack.setChannel(2);
		setMetadata("Label", label_array[a] + " - Angle");
		setPixel(x, y, parseFloat(List.get("Angle")));
		Stack.setChannel(3);
		setMetadata("Label", label_array[a] + " - Area");
		setPixel(x, y, neurite_area);
		close("plot (green)");		
	}
	close("Polar plots");
	close("FFT1");
	selectWindow("Statistics");
	run("Select None");
	selectWindow("canny phi");
	run("Select None");
}