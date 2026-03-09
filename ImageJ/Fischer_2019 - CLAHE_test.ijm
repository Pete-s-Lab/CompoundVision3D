// Macros for applying the CLAHE plugin at http://rsb.info.nih.gov/ij/plugins/clahe/index.html
// by Michael Cammer March 2012

//===================================================================
// This macro operates applies a few different CLAHE paramenters to the same image and pops the
// results into a stack. The paramenters are printed in the upper left of the image. To turn this off, comment out the
// "drawString" line at the end of the loop. On large images, beast run on a subset of the stack. Then pick the
// favorite result and run on the full image.
//===================================================================
setForegroundColor(255, 255, 255);
macro "iterate through CLAHE options" {
original = getImageID;
blocksize = newArray(7, 15, 31, 63, 127, 512);
contrast = newArray(3, 7, 11, 16);
for (k=0; k<contrast.length; k++)
for (i=0; i<blocksize.length; i++) {
selectImage(original);
if (nSlices>1) setSlice(1);
run("Select All");
run("Copy");
run("Add Slice");
run("Paste");
run("Select None");
run("Enhance Local Contrast (CLAHE)", "blocksize="+blocksize[i]+" histogram=256 maximum="+contrast[k]+"");
drawString("blocksize "+blocksize[i]+", contrast "+contrast[k], 1, 30);
}
} // test CLAHE parameters