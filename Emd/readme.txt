Credits go to Yossi Rubner for inventing the emd, which is a wasserstein distance for histogram comparasion, and building the core of this code ;
and Nicolas Lutz, for allowing structures to be allocated on the heap instead of the stack, as well as connecting the dots with ASTex by making the Histogram class.
Original code: http://robotics.stanford.edu/~rubner/emd/default.htm

HOW TO USE:

First, use the ASTex::Histogram class to build two quantized histogram.

To build a quantized histogram:

=============================================================

ImageRGBd im_in;

//Do something with im_in (load, for example)

ImageRGBd::PixelType zero, one; //Specify lowest possible value as zero, and highest as one
for(int i=0; i<3; ++i)
{
	zero[i]=0.0; 
	one[i]=1.0; //could be 255 depending on your type (here, I use 3x double)
}
ASTex::HistogramRGBd histo(im_in);
int nb_bins_per_dimension = 32; //specifies the number of bins per dimension you want to have
ASTex::HistogramRGBBase<int> histo_quantized(zero, one, nb_bins_per_dimension); 

//ASTex::HistogramRGBBase<int> histo_quantized(zero, one); //You can also let the program decide
//The number of bins will then depend on the sizes of the images. 
//The EMD distance works for different quantizations, but is much less precise.
histo_quantized.save("/home/you/some_folder/some_file.csv", nb_bins_per_dimension); //the second argument is new and bound to change, it's necessary for normalizing the wasserstein distance.
//Done

=============================================================

Then, compile the emd program using make, there should really be no problem given its size (and please contact me if there is).
You can then use ./main to compare two histogram files : 

./main example1.csv example2.csv

If you did everything correctly, the result will be a floating point number normalized between 0 and 1, 0 being the result of a comparasion between identical histograms, and 1 being the result of a comparasion between a full black image and a full white image.

Histogram files are as follow :

=============================================================

<Image size, a.k.a. total weight of the bins>
<Number of different bins>
<Number of different bins per dimension (not all of which are necessarily shown in the file)>
<bin_1.x bin_1.y bin_1.z> <weight 1>
..
<bin_n.x bin_n.y bin_n.z> <weight n>
EOF

=============================================================

note that if two histograms compared have too many bins, the measurement can get very slow, yet held the same results as if they were more quantized. 
We suggest maxing at MOST at 32 bins per dimension for RGB images, and depending on how spread the histograms are on the RGB cube, this might be too much (32^3). 
You can see how many different bins your histogram have by looking at the SECOND number in the files.

The program doesn't check yet that what you have fed it is alright, so if you want to make your own histogram, carefully copy the given structure.

Finally, the comparasion doesn't take in account human vision. 
If you want, you can change the float dist(feature_t *F1, feature_t *F2) function in main.c to include, say, deltaE distances.
