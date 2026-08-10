/* Extract regrions of interest (ROIs) for both eyes and the head
 * 
 * This script guides you through the process of extracting one ROI
 * for each eye and one for the whole head.
 * 
 * Peter T. Rühr (October, 2024)
 */

script_version = "0.9.9023";

run("Collect Garbage");
requires("1.39l");
if (isOpen("Log")) { 
     selectWindow("Log"); 
     run("Close"); 
} 
if (isOpen("Results")) { 
     selectWindow("Results"); 
     run("Close"); 
}
if (isOpen("3D Viewer")) { 
     selectWindow("3D Viewer"); 
     run("Close"); 
}
//while (nImages>0) { 
//          selectImage(nImages); 
//          close(); 
//}

plugins = getDirectory("plugins");
unix = '/plugins/';
windows = '\\plugins\\';

if(endsWith(plugins, unix)){
	print("Running on Unix...");
	dir_sep = "/";
}
else if(endsWith(plugins, windows)){
	print("Running on Windows...");
	dir_sep = "\\";
}

//get source dir from user and define other directories
parent_dir_path = getDirectory("Select source Directory");
parent_dir_name = File.getName(parent_dir_path);
specimen_name = parent_dir_name;

log_dir_path = parent_dir_path+dir_sep+'crop_log';
while(File.exists(log_dir_path)){
	waitForUser("Please delete or rename the ROI folder ("+log_dir_path+")!");
}

//filelist = getFileList(parent_dir_path);
// Array.print(filelist);

function listFilesNonRecurse(dir) {
    list = getFileList(dir);
    filelist = newArray();
    for (i = 0; i < list.length; i++) {
        if (!endsWith(list[i], "/")) {
            filelist = Array.concat(filelist, newArray(list[i]));
        }
    }
    return filelist;
}

filelist = listFilesNonRecurse(parent_dir_path);
Array.print(filelist);


if(File.exists(parent_dir_path+dir_sep+"DICOMDIR")){
	print("Trying to open: "+parent_dir_path+"DICOMDIR...");
	run("Bio-Formats", "open=["+parent_dir_path+"DICOMDIR] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); // display_metadata use_virtual_stack
} 
else if (endsWith(filelist[0], "dcm")){
	print("Trying to open dcm files in "+parent_dir_path+dir_sep+"...");
	 run("Image Sequence...", "open="+parent_dir_path+dir_sep+"*.dcm file=.dcm sort");
	File.openSequence(parent_dir_path+dir_sep);
}
else if (endsWith(filelist[0], "tif")){
	if(filelist.length == 1){
		print("Trying to open single tif file in "+parent_dir_path+dir_sep+"...");
		open(parent_dir_path+dir_sep+filelist[0]);
	}
	else {
		print("Trying to open multiple tif files in "+parent_dir_path+dir_sep+"...");
		run("Image Sequence...", "open="+parent_dir_path+dir_sep+"*.tif file=.tif sort");
	}
}
else if (endsWith(filelist[0], "tiff")){
	print("Trying to open multiple tiff files in "+parent_dir_path+dir_sep+"...");
	File.openSequence(parent_dir_path);
}
else if (endsWith(filelist[0], "jp2")){
	print("Trying to open jp2 files in "+parent_dir_path+dir_sep+"...");
	File.openSequence(parent_dir_path);
}
else if (endsWith(filelist[0], ".am")){
	print("Trying to open an Amira(R) file in "+parent_dir_path+dir_sep+"...");
	open(parent_dir_path+dir_sep+filelist[0]);
}
else if (endsWith(filelist[0], ".nrrd")){
	print("Trying to open nrrd file in "+parent_dir_path+dir_sep+"...");
	File.openSequence(parent_dir_path);
	run("Nrrd ...", "load=[+parent_dir_path+dir_sep+filelist[0]]");
}
else if (endsWith(filelist[0], ".czi")){
	print("Trying to open Carl Zeiss Image file in "+parent_dir_path+dir_sep+"...");
	File.openSequence(parent_dir_path);
	open(parent_dir_path+dir_sep+filelist[0]);
} else {
	print("Nothing found in "+parent_dir_path+dir_sep+"...");
}

//OIST:
//run("XRM/TXRM/TXM", "load=[Z:/Jocelyn/Lycaenidae //(adults)/CC_ES003_Arh_bir_F2_Head_IOD/CC_ES003_Arh_bir_F2_Head_IOD_2024-11-20_153623/tomo-//A/CC_ES003_Arh_bir_F2_Head_IOD_tomo-A_recon.txm]");

Stack.getDimensions(width,height,channels,slices,frames);
setSlice(slices/2);
makeRectangle(width/4, height/4, width/2, height/2);
resetMinAndMax();
	
getPixelSize(unit, px_size, ph, pd);

rename("original");
selectWindow("original");
// make sure pixel size is identical in all 3 dimensions
run("Properties...", "channels=1 slices=slices frames=1 pixel_width=px_size pixel_height=px_size voxel_depth=px_size");

// get min and max gray values
run("Set Measurements...", "min redirect=None decimal=0");
run("Measure");
//getResults("Min", 1);
Min=getResult("Min");
Max=getResult("Max");
// print(Min);
// print(Max);
if (isOpen("Results")) { 
     selectWindow("Results"); 
     run("Close"); 
}

// define rectangle for BC
setTool("rectangle");
waitForUser("Define rectangle and get grey values for contrast enhancement.");
run("Brightness/Contrast...");
run("Select None");
Dialog.create("Settings");
Dialog.addMessage("___________________________________");
	Dialog.addNumber("Minimum value: ", Min, 0, 5, " of histogram");
	Dialog.addNumber("Maximum value: ", Max, 0, 5, " of histogram");
	Dialog.show();
	
	min = Dialog.getNumber();
	max = Dialog.getNumber();
	
print("Shifting grey value histogram ...");
//getRawStatistics(count, mean, min, max, std);
run("Select None");
run("Min...", "value=min stack");
run("Max...", "value=max stack");

// define rectangle for head crop
setTool("rectangle");
waitForUser("Define rectangle for head ROI.");
getSelectionCoordinates(x_head, y_head);

waitForUser("1) Check stack for first and last slice number in z direction of head ROI\n2) AFTERWARDS, click 'Ok'."); 
	curr_slice = getSliceNumber();
	Dialog.create("First and last slice");
	Dialog.addMessage("Please enter number of first and last image slice of head ROI.");
	Dialog.addMessage("___________________________________");
	Dialog.addNumber("First image:", 1);
	Dialog.addNumber("Last image:", curr_slice);
	Dialog.show();
	
	first_image_head = Dialog.getNumber();
	last_image_head = Dialog.getNumber();
run("Select None");

// define rectangle for eye1 crop
setTool("rectangle");
waitForUser("Define rectangle for eye1 ROI.");
getSelectionCoordinates(x_eye1, y_eye1);

waitForUser("1) Check stack for first and last slice number in z direction of eye1 ROI\n2) AFTERWARDS, click 'Ok'."); 
	curr_slice = getSliceNumber();
	Dialog.create("First and last slice");
	Dialog.addMessage("Please enter number of first and last slice of eye1 ROI.");
	Dialog.addMessage("___________________________________");
	Dialog.addNumber("First image:", 1);
	Dialog.addNumber("Last image:", curr_slice);
	Dialog.show();
	
	first_image_eye1 = Dialog.getNumber();
	last_image_eye1 = Dialog.getNumber();

run("Select None");

// define rectangle for eye2 crop
setTool("rectangle");
waitForUser("Define rectangle for eye2 ROI.");
getSelectionCoordinates(x_eye2, y_eye2);

waitForUser("1) Check stack for first and last slice number in z direction of eye2 ROI\n2) AFTERWARDS, click 'Ok'."); 
	curr_slice = getSliceNumber();
	Dialog.create("First and last slice");
	Dialog.addMessage("Please enter number of first and last slice id eye2 ROI.");
	Dialog.addMessage("___________________________________");
	Dialog.addNumber("First image:", 1);
	Dialog.addNumber("Last image:", curr_slice);
	Dialog.show();
	
	first_image_eye2 = Dialog.getNumber();
	last_image_eye2 = Dialog.getNumber();


// check if scaling is necessary for head-stack
Dialog.create("Scaling settings");
	Dialog.addMessage("___________________________________");
	Dialog.addNumber("Scale to [MB]: ", 120);
	Dialog.addMessage("___________________________________");
	Dialog.addMessage("PTR, Jan. 2022-2024");
	Dialog.addMessage("Bonn Univ., Bonn, Germany");
	Dialog.show();
	d_size = Dialog.getNumber()/1024;  //MB/1024=GB

//Create dialog to check if pixel size is correct or user define it
Dialog.create("Check pixel size");
	Dialog.addNumber("Correct pixel size?:", px_size, 9, 15, unit)
	Dialog.show();
	px_size = Dialog.getNumber();
	unit = Dialog.getString();

// set pixil size
run("Properties...", "pixel_width="+px_size+" pixel_height="+px_size+" voxel_depth="+px_size);

// deselct everything and convert stack to 8bit
run("Select None");
run("8-bit");

// crop head
selectWindow("original");
run("Duplicate...", "duplicate range=first_image_head-last_image_head");
makeRectangle(x_head[0], y_head[0], x_head[2]-x_head[0], y_head[2]-y_head[0]);
run("Crop");
title_head = getTitle();

//Create target directory to save tiffs/nrrds in
print("Creating target directory...");
ROI_name = "head";
ROI_path = parent_dir_path+dir_sep+specimen_name+"_"+ROI_name;
print(ROI_path);
File.makeDirectory(ROI_path);
// sace full-sized head stack
// saveAs("Tiff", ROI_path+dir_sep+specimen_name+"_"+ROI_name+".tif");

// calculate if scaling is necessary later
Stack.getDimensions(width_orig, height_orig, channels, slices, frames);
o_size = width_orig*height_orig*slices/(1024*1024*1024);
print("Target directory loaded has a size of ~"+o_size+" GB.");
d = pow(d_size/o_size,1/3);
perc_d = round(100 * d);
d = perc_d/100;

if(perc_d<100){
	print("Scaling head stack to "+perc_d+"% to reach stack size of ~"+d_size+" GB...");
	run("Scale...", "x="+d+" y="+d+" z="+d+" interpolation=Bicubic average process create");
	print("New px size = "+px_size+" um.");
	// saveAs("Tiff", ROI_path+dir_sep+specimen_name+"_"+ROI_name+"_sc"+perc_d+".tif");
	nrrd_file = ROI_path+dir_sep+specimen_name+"_"+ROI_name+"_sc"+perc_d+".nrrd";
	run("Nrrd ... ", "nrrd=[nrrd_file]");
} else {
	print("No scaling of head stack necessary.");
	perc_d = 100;
	nrrd_file = ROI_path+dir_sep+specimen_name+"_"+ROI_name+".nrrd";
	run("Nrrd ... ", "nrrd=[nrrd_file]");
}

getPixelSize(unit_, px_size_sc, ph, pd);


// crop eye 1
selectWindow("original");
run("Duplicate...", "duplicate range=first_image_eye1-last_image_eye1");
//run("Make Substack...", " slices="+first_image_eye1+"-"+last_image_eye1);
makeRectangle(x_eye1[0], y_eye1[0], x_eye1[2]-x_eye1[0], y_eye1[2]-y_eye1[0]);
run("Crop");
title_eye1 = getTitle();
print(title_eye1 + " -> " + "eye1");
rename("eye1");
title_eye1 = getTitle();

//Create target directory to save tiffs/nrrds in
print("Creating target directory...");
ROI_name = "eye1";
ROI_path = parent_dir_path+dir_sep+specimen_name+"_"+ROI_name;
File.makeDirectory(ROI_path);
//run("Image Sequence... ", "dir=" + ROI_path + " format=TIFF name=" + specimen_name+"_"+ROI_name+"_");
// saveAs("Tiff", ROI_path+dir_sep+specimen_name+"_"+ROI_name+".tif");
nrrd_file = ROI_path+dir_sep+specimen_name+"_"+ROI_name+".nrrd";
run("Nrrd ... ", "nrrd=[nrrd_file]");

// save STL
STL_path = parent_dir_path+dir_sep+specimen_name+"_"+ROI_name+dir_sep+"stl";
File.makeDirectory(STL_path);

Stack.getDimensions(width_eye1,height_eye1,channels_eye1,slices_eye1,frames_eye1);
setSlice(slices_eye1/2);
setTool("line");


// crop eye 2
selectWindow("original");
run("Duplicate...", "duplicate range=first_image_eye2-last_image_eye2");
// run("Make Substack...", " slices="+first_image_eye2+"-"+last_image_eye2);
makeRectangle(x_eye2[0], y_eye2[0], x_eye2[2]-x_eye2[0], y_eye2[2]-y_eye2[0]);
run("Crop");
title_eye2 = getTitle();
print(title_eye2 + " -> " + "eye2");
rename("eye2");
title_eye2 = getTitle();

//Create target directory to save tiffs/nrrds in
print("Creating target directory...");
ROI_name = "eye2";
ROI_path = parent_dir_path+dir_sep+specimen_name+"_"+ROI_name;
File.makeDirectory(ROI_path);
//run("Image Sequence... ", "dir=" + ROI_path + " format=TIFF name=" + specimen_name+"_"+ROI_name+"_");
//saveAs("Tiff", ROI_path+dir_sep+specimen_name+"_"+ROI_name+".tif");
nrrd_file =ROI_path+dir_sep+specimen_name+"_"+ROI_name+".nrrd";
run("Nrrd ... ", "nrrd=[nrrd_file]");

// save STL
STL_path = parent_dir_path+dir_sep+specimen_name+"_"+ROI_name+dir_sep+"stl";
File.makeDirectory(STL_path);

Stack.getDimensions(width_eye2,height_eye2,channels_eye2,slices_eye2,frames_eye2);
setSlice(slices_eye2/2);

Dialog.create("Extract surface meshes?");
Dialog.addMessage("___________________________________");
	Dialog.addCheckbox("Extract surface meshes", true);
	Dialog.show();
	
	extract_surfaces = Dialog.getCheckbox();
	
if(extract_surfaces == true){
	selectWindow(title_eye1);
	waitForUser("Threshold", "Find appropriate threshold level for eye1.");
		Dialog.create("Threshold");
		Dialog.addMessage("Please set appropriate threshold level for eye1.");
		Dialog.addNumber("Threshold :", 100);
		Dialog.show();
		
		threshold_eye1 = Dialog.getNumber();
	
	
	selectWindow(title_eye2);
	waitForUser("Threshold", "Find appropriate threshold level for eye2.");
		Dialog.create("Threshold");
		Dialog.addMessage("Please set appropriate threshold level for eye2.");
		Dialog.addNumber("Threshold :", threshold_eye1);
		Dialog.show();
		
		threshold_eye2 = Dialog.getNumber();
} else {
	threshold_eye1 = "NA";
	threshold_eye2 = "NA";
}


print("************************************");
File.makeDirectory(log_dir_path);
ROI_log_file = specimen_name+"_crop.log";
print("Creating log file ("+ROI_log_file+")...");
open_log_file = File.open(log_dir_path+dir_sep+ROI_log_file);
print(open_log_file, "script_version = " + script_version);
print(open_log_file, "px_size = " + px_size + " " + unit);
print(open_log_file, "hist_min = " + min);
print(open_log_file, "hist_max = " + max);
print(open_log_file, "ROI_head = makeRectangle("+x_head[0]+", "+y_head[0]+", "+x_head[1]-x_head[0]+", "+y_head[2]-y_head[0]+");");
print(open_log_file, "z_first_head = " + first_image_head);
print(open_log_file, "z_last_head = " + last_image_head);
print(open_log_file, "ROI_eye1 = makeRectangle("+x_eye1[0]+", "+y_eye1[0]+", "+x_eye1[1]-x_eye1[0]+", "+y_eye1[2]-y_eye1[0]+");");
print(open_log_file, "z_first_eye1 = " + first_image_eye1);
print(open_log_file, "z_last_eye1 = " + last_image_eye1);
print(open_log_file, "threshold_eye1 = " + threshold_eye1);
print(open_log_file, "ROI_eye2 = makeRectangle("+x_eye2[0]+", "+y_eye2[0]+", "+x_eye2[1]-x_eye2[0]+", "+y_eye2[2]-y_eye2[0]+");");
print(open_log_file, "z_first_eye2 = " + first_image_eye2);
print(open_log_file, "z_last_eye2 = " + last_image_eye2);
print(open_log_file, "threshold_eye2 = " + threshold_eye2);
print(open_log_file, "scaling_head = " + d);
print(open_log_file, "px_size_head = " + px_size_sc + " " + unit);
File.close(open_log_file)
print("************************************");

if(extract_surfaces == true){
	print("Extracting surfaces");
	selectWindow(title_eye1);
	run("3D Viewer");
	call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
	call("ij3d.ImageJ3DViewer.add", title_eye1, "White", title_eye1, threshold_eye1, "true", "true", "true", "1", "2");
	
	waitForUser("Save STL...", "Save binary STL file for eye1 manually. Click okay AFTER saving is finished.");
	call("ij3d.ImageJ3DViewer.close");
	
	
	selectWindow(title_eye2);
	run("3D Viewer");
	call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
	call("ij3d.ImageJ3DViewer.add", title_eye2, "White", title_eye2, threshold_eye2, "true", "true", "true", "1", "2");
	
	waitForUser("Save STL...", "Save binary STL file for eye2 manually. Click okay AFTER saving is finished.");
	call("ij3d.ImageJ3DViewer.close");
}

while (nImages>0) { 
          selectImage(nImages); 
          close(); 
}

run("Collect Garbage");
print("All done!");
print("************************************");