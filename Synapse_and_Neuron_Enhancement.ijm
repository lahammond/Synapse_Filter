// CIZI 3D Cell Counter
//
// Author: 	Luke Hammond
// Cellular Imaging | Zuckerman Institute, Columbia University
// Version: 0.1
// Date:	4th April 2018
//
// For detection of cells in 3D images. Generates a CSV file containing total counts and image #
// Expects multichannel TIF from 3D Brain Region Enhancement


// Initialization
requires("1.51w");



#@ File[] listOfPaths(label="select files or folders", style="both")

#@ boolean(label="Tiled datasets: stitch tiles after processing?", description="Enable for stitching post processing.") StitchON
#@ Integer(label="Grid size X:", value = 3, style="spinner", description="") XTiles
#@ Integer(label="Grid size X:", value = 3, style="spinner", description="") YTiles
#@ Integer(label="Resolution used for processing (0 for original resolution):", value = 0, style="spinner", description="") FinalRes

#@ String(label="Select syanpse channel for filtering:", choices={"1", "2", "3", "4", "0"}, style="radioButtonHorizontal", value = "1", description="") CellCh
#@ Integer(label="Unsharp Mask radius (px, 0 if none):", value = 1, style="spinner", description="") USSZ
#@ BigDecimal(label="Unsharp Mask weight:", value = 0.700, style="spinner") USMW
#@ Integer(label="Background removal filter radius (px):", value = 5, style="spinner", description="") OSBSRad
#@ Integer(label="Background removal filter iterations:", value = 20, style="spinner", description="") OSBSItr
// Integer(label="Post filter remove outliers:", value = 4, style="spinner", description="") RORad
#@ boolean(label="3D median filter", description="1px median filter recommened for 3D stacks.") Med3DON


#@ String(label="Select neuron channel for filtering (0 if none):", choices={"1", "2", "3", "4", "0"}, style="radioButtonHorizontal", value = "0", description="") CellCh2
#@ boolean(label="3D median filter", description="1px median filter recommened for 3D stacks.") Med3DON2



starttime = getTime();
run("Options...", "iterations=3 count=1 black do=Nothing");
run("Set Measurements...", "fit redirect=None decimal=3");
run("Colors...", "foreground=white background=black selection=yellow");
run("Clear Results");


setBatchMode(true);


// Preparation


for (FolderNum=0; FolderNum<listOfPaths.length; FolderNum++) {
	input=listOfPaths[FolderNum];
	if (File.exists(input)) {
        if (File.isDirectory(input)) {
			run("Clear Results");
			print("\\Clear");
			print("Synapse and neuron enhancement running:");
			print(" Processing folder "+(FolderNum+1)+" of "+listOfPaths.length+" folders selected for processing." );
	        print(" Processing folder "+(FolderNum+1)+": " + input + " ");

			// Create folders
			File.mkdir(input + "/Processed");
			ChOut = input + "/Processed/";
			Ch2Out = input + "/Processed/";
			// Create folders
			File.mkdir(input + "/Processed_MIPs");
			ChMIPOut = input + "/Processed_MIPs/";
			Ch2MIPOut = input + "/Processed_MIPs/";

/*
File.mkdir(input + "Analyzed/Channel"+CellCh+"_Object_Validation");
ChOut = input + "Analyzed/Channel"+CellCh+"_Object_Validation/";

if (CellCh2 > 0) {
	File.mkdir(input + "Analyzed/Channel"+CellCh2+"_Object_Validation");
	Ch2Out = input + "Analyzed/Channel"+CellCh2+"_Object_Validation/";
}


//Create Table - either have open the whole time or open and close after each iteration
TableTitle = "Cell_Counts"; 
TableTitle = "["+TableTitle+"]"; 
f=TableTitle;  
run("New... ", "name="+TableTitle+" type=Table"); 
//CellCh2=2;
//CellCh=1;
if (CellCh2 > 0) {
	print(f,"\\Headings:Image\tCountCh"+CellCh+"\tCountCh"+CellCh2);
} else {
	print(f,"\\Headings:Image\tCountCh"+CellCh);
}

*/
// get files
files = getFileList(input);	
files = ImageFilesOnlyArray(files);		

run("Collect Garbage");

//iterate over all files

for(i=0; i<files.length; i++) {				
	image = files[i];	
	print("\\Update3: Processing image " + (i+1) +" of " + files.length +".");
	run("Bio-Formats Importer", "open=[" + input + "/" + image + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	getPixelSize(unit, W, H);
	getDimensions(width, height, ChNum, slices, frames);
	
	
	//Rescale Image 
	if (FinalRes > 0) {
		RescaleImage();
	} 

	

	rename("Raw");

	// Process First Cell Channel

	// multichannel?
	if (ChNum > 1) {
		run("Split Channels");
		selectWindow("C" + CellCh + "-Raw");
	}
	rename("Ch1");
	
	//run("Duplicate...", "title=Ch1-Raw duplicate");

	selectWindow("Ch1");
	//Process Image
	run("Unsharp Mask...", "radius="+USSZ+" mask="+USMW+" stack");
	run("Unsharp Mask...", "radius="+USSZ+" mask="+USMW+" stack");
	OSBSFilter("Ch1", OSBSRad, OSBSItr);

	//if (RORad > 0) {
	//	run("Remove Outliers...", "radius="+RORad+" threshold=0 which=Bright stack");
	//}


	if (Med3DON == true) {
		run("Median 3D...", "x=1 y=1 z=1");
	}
	
	
	//run("Remove Outliers...", "radius="+RORad+" threshold=0 which=Bright stack");
	//run("Remove Outliers...", "radius=4 threshold=10 which=Dark stack");

	if (CellCh2 == 0 ) {
		saveAs("Tiff", ChOut + "C"+CellCh+"_"+ image);
		rename("Ch1");
		run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", ChMIPOut + "C"+CellCh+"_"+ image);
		close();
		selectWindow("Ch1");
		close();
	} else {
		rename("C1");
	}
	

	
	if (CellCh2 > 0) {
		selectWindow("C" + CellCh2 + "-Raw");
		rename("Ch2");
		
	
		if (Med3DON2 == true) {
			run("Median 3D...", "x=1 y=1 z=1");
		}
	
		rename("C2");
		
		//run("Remove Outliers...", "radius="+RORad2+" threshold=0 which=Bright stack");
		//run("Remove Outliers...", "radius=4 threshold=10 which=Dark stack");
		//saveAs("Tiff", Ch2Out + "C"+CellCh2+"_"+ image);
		
		run("Merge Channels...", "c1=C1 c2=C2 create");
		Stack.setChannel(1);
		run("Green");
		Stack.setChannel(2);
		run("Red");
		if (StitchON == true) {
			if (i <= 8) {
				saveAs("Tiff", Ch2Out + "tile_0"+(i+1)+".tif");
			} else {
				saveAs("Tiff", Ch2Out + "tile_"+(i+1)+".tif");
			}
		} else {
			saveAs("Tiff", Ch2Out + image);	
		}
		
		rename("Merge");
		
		run("Z Project...", "projection=[Max Intensity]");
		if (StitchON == true) {
			if (i <= 8) {
				saveAs("Tiff", Ch2MIPOut + "tile_0"+(i+1)+".tif");
			} else {
				saveAs("Tiff", Ch2MIPOut + "tile_"+(i+1)+".tif");
			}
		} else {
			saveAs("Tiff", Ch2MIPOut + image);
		}
		
		// turned off image name and now saving tiles
		
		close();
		selectWindow("Merge");
		close();

		
		run("Collect Garbage");

		}
		}

		// Stitch dataset
		if (StitchON == true) {	
			run("Grid/Collection stitching", "type=[Grid: snake by rows] order=[Right & Down                ] grid_size_x="+XTiles+" grid_size_y="+YTiles+" tile_overlap=10 first_file_index_i=1 directory=["+Ch2Out+"] file_names=tile_{ii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap subpixel_accuracy display_fusion computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
			saveAs("Tiff", input + "/Stitched_Image.tif");
		}
		// Timing information


        
        endtime = getTime();
		dif = (endtime-starttime)/1000;
		print("  ");
		print("----------------------------------------------------------------");
		print("Synapse and neuron enhancement completed. Processing time =", (dif/60), "minutes. ", (dif/i), "seconds per image.");
		print("----------------------------------------------------------------");
		selectWindow("Log");
		saveAs("txt", input+"/OSBS_Filter_Log.txt");
        }
	}
}







//run("Remove Outliers...", "radius=2 threshold=0 which=Bright stack");
//run("Unsharp Mask...", "radius=4 mask=0.70 stack");
//run("Enhance Contrast...", "saturated=0.3 process_all use");
//run("Remove Outliers...", "radius=4 threshold=0 which=Bright stack");








function OSBSFilter(imagename, radius, iterations) {
	selectWindow(imagename);
	getDimensions(w2, h2, c2, slices, f2);
	rename("ROFimage");
	run("Duplicate...", "title=bgstack duplicate");
	run("Z Project...", "projection=[Max Intensity]");
	for(k=0; k<iterations; k++) {	
		run("Remove Outliers...", "radius="+radius+" threshold=0 which=Bright");
		//run("Morphological Filters", "operation=Dilation element=Disk radius=3");
	}
	run("Morphological Filters", "operation=Dilation element=Disk radius=10");
	rename("morph");
	selectWindow("MAX_bgstack");
	close();
	selectWindow("morph");
	rename("MAX_bgstack");
	
	for (m=0; m<slices; m++) {
		selectWindow("MAX_bgstack");
		run("Select All");
		run("Copy");
		selectWindow("bgstack");
		setSlice(m+1);
		run("Paste");
	}
	selectWindow("MAX_bgstack");
	close();
	
	imageCalculator("Subtract create stack", "ROFimage","bgstack");
	selectWindow("bgstack");
	close();
	selectWindow("ROFimage");
	close();
	selectWindow("Result of ROFimage");
	rename(imagename);
}	

function ImageFilesOnlyArray (arr) {
	//pass array from getFileList through this e.g. NEWARRAY = ImageFilesOnlyArray(NEWARRAY);
	setOption("ExpandableArrays", true);
	f=0;
	files = newArray;
	for (i = 0; i < arr.length; i++) {
		if(endsWith(arr[i], ".tif") || endsWith(arr[i], ".nd2") ) {   //if it's a tiff image add it to the new array
			files[f] = arr[i];
			f = f+1;
		}
	}
	arr = files;
	arr = Array.sort(arr);
	return arr;
}

function DeleteDir(Dir){
	listDir = getFileList(Dir);
  	//for (j=0; j<listDir.length; j++)
      //print(listDir[j]+": "+File.length(myDir+list[i])+"  "+File. dateLastModified(myDir+list[i]));
 // Delete the files and the directory
	for (j=0; j<listDir.length; j++)
		ok = File.delete(Dir+listDir[j]);
	ok = File.delete(Dir);
	if (File.exists(Dir))
	    print("\\Update10: Unable to delete temporary directory"+ Dir +".");
	else
	    print("\\Update10: Temporary directory "+ Dir +" and files successfully deleted.");
}

function RescaleImage(){
	//Expects FinalRes as an input from user in menu
	input_Title = getTitle();
	input_ID = getImageID();
	//get image information		
	getPixelSize(unit, W, H);
	// Determine rescale value
	Rescale = (1/(FinalRes/W));
	run("Scale...", "x="+Rescale+" y="+Rescale+" interpolation=Bilinear average create");
	rescale_ID = getImageID(); 
	selectImage(input_ID);
	close();
	selectImage(rescale_ID);
	rename(input_Title);
}

function NumberedArray(maxnum) {
	//use to create a numbered array from 1 to maxnum, returns numarr
	//e.g. ChArray = NumberedArray(ChNum);
	numarr = newArray(maxnum);
	for (i=0; i<numarr.length; i++){
		numarr[i] = (i+1);
	}
	return numarr;
}

function closewindow(windowname) {
	if (isOpen(windowname)) { 
      		 selectWindow(windowname); 
       		run("Close"); 
  		} 
}