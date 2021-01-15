close("*");

//Filter variables
axial_blur = 2; //Size of axial Gaussian blur
lateral_median = 3; //Size of lateral median filter to suppress neurites
gamma = 0.6; //Gamma applied to final image

//Direction variables
//FFT_substack_nSlices = 3; //Thickness in slices of substack to generate FFT (ensures neurites on different slices are included);
//n_stdev = 2; //Number of standard deviations above the mean to threshold the FFT
tolerance_range = 0.1; //+/- factor range to test tolerance
tolerance_step = tolerance_range; //Step size to test tolerance
min_size_cutoff = 4; //Smallest size to be exluded in particle analyzer
max_size_cutoff = 20; //Largest size to be exluded in particle analyzer
FFT_min_radius = 1/16; //Ratio of total width of FFT to set as min cutoff - i.e. lowest sptatial frequency
FFT_max_radius = 1/4; //Ratio of total width of FFT to set as max cutoff - i.e. highest sptatial frequency
sector_width = 10; //Width of sector off FFT to look for directionality
n_degrees = 360; //Number of degrees to sample in FFT

//Initialize global variables
tolerance = -1;
size = -1;
label_array = newArray("Mean", "Median", "Min", "Max");

//Process image to enrichj for 
setBatchMode(true);
process_files = true;
out_dir = getDirectory("Choose an output directory:");
file_list = getFileList(out_dir);

while(process_files){
	close("*");
	file_path = File.openDialog("Select file to process:");
	for(a=0; a<file_list.length; a++){ //Check if file has already been processed
		if(startsWith(file_list[a], "Statistics - ")){
			sample = replace(file_list[a], "Statistics - ", "");
			sample = replace(sample, ".tif", "");
			if(matches(file_path, ".*" + sample + ".*")) process_files = getBoolean("Sample: \"" + sample + "\" has already been processed, do you want to process it again?");
		}
	}
	if(process_files){
		run("Bio-Formats Importer", "open=[" + file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		i = getTitle();
		newImage("Statistics", "32-bit black", n_degrees, (2*label_array.length)+1, 1);
		processImage(i);
		traceDendrites("blur");
		getDirectionality(i, tolerance);
	}
	process_files = getBoolean("Would you like to process another file?");
}
updateResults();
if(isOpen("Results")) {
	selectWindow("Results");
	run("Close");
}
close("*");
setBatchMode("exit and display");
exit();

function processImage(i){
	selectWindow(i);
	run("Split Channels");
	selectWindow("C2-" + i);
	run("Duplicate...", "title=blur duplicate");
	selectWindow("blur");
	run("Gaussian Blur 3D...", "x=0 y=0 z=" + axial_blur); //Axial Gaussian blur to remove high temporal frequency noise
	selectWindow("blur");
	run("Duplicate...", "title=median duplicate");
	selectWindow("median");
	run("Median...", "radius=" + lateral_median + " stack"); //Lateral median to get background enriched image
	imageCalculator("Subtract stack", "blur","median"); //Subtract median from previous image to get image enriched for neurites
	close("median");
	selectWindow("blur");
	run("Gamma...", "value=" + gamma + " stack"); //Apply gamma to make dim peojections easier to see
	selectWindow("blur");
	close("C1-" + i);
	close("C3-" + i);
}

//This function traces dendrites based on the local maxima of the XZ and YZ planes
function traceDendrites(image){
	accept_tolerance = false;
	tolerance = -1;
	size_cutoff = -1;
	while(!accept_tolerance){
		//Find tolernace threshold for dendrites
		selectWindow(image);
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("MAX_blur");
		setBatchMode("show");
		run("Threshold...");
		if(tolerance == -1) setAutoThreshold("Triangle dark");
		else setThreshold(tolerance, upper);
		waitForUser("Please set the threshold and then click \"OK\"");
		getThreshold(tolerance, upper);
		selectWindow("Threshold");
		run("Close");
		close("MAX_blur");
		selectWindow(image);
		
		//Store the image voxel dimensions
		getVoxelSize(voxelWidth, voxelHeight, voxelDepth, voxelUnit);
	
		//Temporarily convert voxels to pixel dimensions to prevent pixel interpolation during transform
		setVoxelSize(1, 1, 1, "pixel");

//-----------First, trace the dendrites in the YZ plane--------------------
		concat_string = "  title=[Threshold stack] open ";
		n_stack = 0;
		for(t=-1*tolerance_range; t<=tolerance_range; t+=tolerance_step){
			//Rotate the image to get a YZ projection
			n_stack++;
			selectWindow(image);
			run("TransformJ Rotate", "z-angle=0.0 y-angle=270 x-angle=0.0 interpolation=linear background=0.0 adjust");
			run("Rotate 90 Degrees Left");
		
			//Find local maxima in each slice (dendrites)
			selectWindow(image + " rotated");
			stackLength = nSlices;
		
			for(a=1; a<=stackLength; a++){
				selectWindow(image + " rotated");
				setSlice(a);
				run("Find Maxima...", "noise=" + round(tolerance*(1+t)) + " output=[Single Points]");
				selectWindow(image + " rotated(" + a + ") Maxima");
				if(a==1) rename("YZ_stack");
				else run("Concatenate...", "  title=[YZ_stack] image1=[YZ_stack] image2=[" + image + " rotated(" + a + ") Maxima] image3=[-- None --]");
			}
		
			//Close the original image so the points can be concatenated
			close(image + " rotated");
		
			//Rotate the stack back to its original orientation
			selectWindow("YZ_stack");
			run("Rotate 90 Degrees Right");
			run("TransformJ Rotate", "z-angle=0.0 y-angle=90 x-angle=0.0 interpolation=linear background=0.0 adjust");
			close("YZ_stack");
			selectWindow("YZ_stack rotated");
			rename("YZ_stack");
		
			//Rotate the image to get an XZ projection
			selectWindow(image);
			run("TransformJ Rotate", "z-angle=0.0 y-angle=0.0 x-angle=270 interpolation=linear background=0.0 adjust");
			run("Rotate 90 Degrees Left");
			run("Rotate 90 Degrees Left");
		
			//Find local maxima in each slice (dendrites)
			selectWindow(image + " rotated");
			stackLength = nSlices;
			for(a=1; a<=stackLength; a++){
				selectWindow(image + " rotated");
				setSlice(a);
				run("Find Maxima...", "noise=" + round(tolerance*(1+t)) + " output=[Single Points]");
				selectWindow(image + " rotated(" + a + ") Maxima");
				if(a==1) rename("XZ_stack");
				else run("Concatenate...", "  title=[XZ_stack] image1=[XZ_stack] image2=[" + image + " rotated(" + a + ") Maxima] image3=[-- None --]");
			}
		
			//Close the original image so the points can be concatenated
			close(image + " rotated");
		
			//Rotate the stack back to its original orientation
			selectWindow("XZ_stack");
			run("Rotate 90 Degrees Right");
			run("Rotate 90 Degrees Right");
			run("TransformJ Rotate", "z-angle=0.0 y-angle=0.0 x-angle=90 interpolation=linear background=0.0 adjust");
			close("XZ_stack");
			selectWindow("XZ_stack rotated");
			rename("XZ_stack");
		
			//Add the XZ and YZ dendrite maps to get the whole map
			imageCalculator("Add stack", "YZ_stack","XZ_stack");
			selectWindow("YZ_stack");
			rename("Dendrite mask");
			close("XZ_stack");
		
			//Remove disconnected points
			selectWindow("Dendrite mask");
			run("Z Project...", "projection=[Max Intensity]");
			for(a=min_size_cutoff; a<=max_size_cutoff; a++){
				selectWindow("MAX_Dendrite mask");
				run("Analyze Particles...", "size=" + a + "-Infinity show=Masks");
				selectWindow("Mask of MAX_Dendrite mask");
				if(a == min_size_cutoff) rename("Particle Analyzer Stack - tolerance " + round(tolerance*(1+t)));
				else{
					run("Concatenate...", "  title=[Particle Analyzer Stack - tolerance " + round(tolerance*(1+t)) + "] image1=[Particle Analyzer Stack - tolerance " + round(tolerance*(1+t)) + "] image2=[Mask of MAX_Dendrite mask] image3=[-- None --]");
				}
			}
			selectWindow("Dendrite mask");
			rename("Dendrite mask" + n_stack);
			close("Max_Dendrite mask");
			concat_string += "image" + n_stack + "=[Particle Analyzer Stack - tolerance " + round(tolerance*(1+t)) + "] "; 
		}
		//Create hyperstack of toelrances and sizes to allow user to select optimal setting
		if(n_stack > 1) run("Concatenate...", concat_string);
		else{
			selectWindow("Particle Analyzer Stack - tolerance " + tolerance);
			rename("Threshold stack");
		}
		selectWindow("Threshold stack");
		Stack.getDimensions(width, height, channels, slices, frames);
		newImage("Blur stack", "8-bit color-mode", width, height, 1, slices, frames);
		selectWindow("blur");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("MAX_blur");
		run("8-bit");
		run("Select All");
		run("Copy");
		selectWindow("Blur stack");
		for(f=1; f<=frames; f++){
			for(s=1; s<=slices; s++){
				Stack.setSlice(s);
				Stack.setFrame(f);
				run("Paste");
			}
		}
		run("Select None");
		run("Merge Channels...", "c2=[Blur stack] c6=[Threshold stack] create ignore");
		if(isOpen("Composite")){ //1 frame stacks are called composite on merge - hyperstacks are called merge
			selectWindow("Composite");
			rename("Merged");
		}
		selectWindow("Merged");
		setBatchMode("show");
		waitForUser("Select the optimal size (Z) and threshold (t)");
		

		//Create overlay of mask and image to check if tolerance is valid
		setBatchMode("hide");
		selectWindow("Merged");
		Stack.getPosition(dummy, size_slice, tolerance_frame);
		tolerance = round(tolerance * ((1-tolerance_range) + ((tolerance_frame-1) * tolerance_step)));
		size = size_slice + min_size_cutoff - 1;
		run("Duplicate...", "title=[Mask of MAX_Dendrite mask] duplicate channels=2 slices=size_slice frames=tolerance_frame");
		selectWindow("Merged");
		run("Duplicate...", "title=[Mask of MAX_Dendrite mask2] duplicate channels=2 slices=size_slice frames=tolerance_frame");
		close("Merged");
		run("Merge Channels...", "c2=MAX_blur c6=[Mask of MAX_Dendrite mask2] create ignore");
		selectWindow("Composite");
		setBatchMode("show");
		for(a=0; a<4; a++){
			Stack.setDisplayMode("color");
			wait(1000);
			Stack.setDisplayMode("composite");
			wait(1000);
		}
		accept_tolerance = getBoolean("Accept tolerance threshold? (Tolerance = " + tolerance + ", size = " + size + ")");
		if(!accept_tolerance){
			for(f=1; f<=frames; f++) close("Dendrite mask" + f);
		}
		else{
			imageCalculator("AND stack", "Dendrite mask" + tolerance_frame,"Mask of MAX_Dendrite mask");
			for(f=1; f<=frames; f++){
				selectWindow("Dendrite mask" + f);
				if(f == tolerance_frame) rename("Dendrite mask");
				else close("Dendrite mask" + f);
			}
		}
		close("Composite");
	}
	//Add size and tolerance to output matrix
	selectWindow("Statistics");
	setPixel(0, 2*label_array.length, tolerance); 
	setPixel(1, 2*label_array.length, size); 
}

function getDirectionality(image, tolerance){
	//Get FFT of mask
	run("FFT Options...", "raw");
	selectWindow("Dendrite mask");
	run("Z Project...", "projection=[Max Intensity]");
	selectWindow("MAX_Dendrite mask");
	run("FFT");
	selectWindow("PS of MAX_Dendrite mask");
	rename("FFT1");
	//run("Log");
	//run("Gaussian Blur...", "sigma=8");

	//make angle map
	selectWindow("FFT1");
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
	close("canny phi");

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
		selectWindow("Statistics");
		setPixel(2, label_array.length + a, ellipse_ratio);
		setPixel(3, label_array.length + a, parseFloat(List.get("Angle")));
		close("plot (green)");		
	}
	selectWindow("Polar plots");
	run("Select None");
	image = replace(image, "\\..+?$", ""); 
	saveAs("tiff", out_dir + "Polar plots - " + image);
	selectWindow("Statistics");
	run("Select None");
	saveAs("tiff", out_dir + "Statistics - " + image);
	selectWindow("Dendrite mask");
	run("Select None");
	saveAs("tiff", out_dir + "Dendrite mask - " + image);
	selectWindow("FFT1");
	run("Select None");
	saveAs("tiff", out_dir + "FFT - " + image);
}






