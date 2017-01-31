# Installation

* Install [StereoPipeline](https://ti.arc.nasa.gov/tech/asr/intelligent-robotics/ngt/stereo/) and add its /bin directory to $PATH

* Install [libLAS](http://www.liblas.org/) : ``sudo apt-get install liblas-bin``

* Compile source with compile.sh

* (optional) For visualization, you can add [CImg.h](http://cimg.eu/) library and add "-DUSE_VIS -lX11" to compilation options.

---

# Command-line Arguments

## Modes

``-v INPUT [OUTPUT] [-vscale SCALE]`` 

Visualization mode. Program has to be compiled with visualization support in order for it to work. Reads point cloud ``INPUT`` file and saves it as ``OUTPUT`` image (only bmp format is supported). If ``OUTPUT`` is omitted then the image will have the same name as ``INPUT`` with adjusted extension. You can add optional ``-vscale`` option to shrink/enlarge color space. Default is 1.0.


``-m 1 INPUT_1 INPUT_2 ... INPUT_N OUTPUT [-lima LIMIT_A] [-limb LIMIT_B] [-limc LIMIT_C] [-check CHECK] [-maug MERGE_AUGMENT_DIST] [-opt] [-exp] [-fill] [-t THREADS_NO]`` 

Merge mode. Reads all of the ``INPUT_X`` point clouds, merges them and saves the merged point cloud to ``OUTPUT`` file. Other optional arguments controls how the merging is done. ``-exp`` adds exploiting of scoring function (used during the contest), ``-fill`` fills small holes in the input point clouds (difference is negligible), ``-opt`` improves speed of the second iteration of merging algorithm (no point in omitting it), ``-t THREADS_NO`` sets number of threads to ``THREADS_NO`` (default is 4, anything above 8 won't improve speed), ``-maug MERGE_AUGMENT_DIST`` adds naive augmentation of each point cloud by creating four additional copies at (-MERGE_AUGMENT_DIST, 0), (+MERGE_AUGMENT_DIST, 0), (0, -MERGE_AUGMENT_DIST), (0, +MERGE_AUGMENT_DIST) offsets and finally ``-lima LIMIT_A``, ``-limb LIMIT_B``, ``-limc LIMIT_C`` and ``-check CHECK`` are various constants used in scoring function exploting.


``-r 1 OUTPUT_PREFIX INPUT_A INPUT_B "+OPTION_LIST" -nitfdir NITF_DIRECTORY -kml KML -data SIMPLIFIED_DATA_DIRECTORY``

Run Mode. Reads ``NITF_DIRECTORY/INPUT_A`` and ``NITF_DIRECTORY/INPUT_B`` ntf files and produces ``results/OUTPUT_PREFIX.txt`` point cloud file. ``KML`` points out to kml file, and ``SIMPLIFIED_DATA_DIRECTORY`` points out to directory as specified in the original problem statement (parent directory that contains rpc, metadata and cropped images). ``OPTION_LIST`` is a space-delimited list of options that are directly pushed to StereoPipeline software. 


``-f 1 OUTPUT_PREFIX ID_1 ID_2 ... ID_N "+OPTION_LIST" -nitfdir NITF_DIRECTORY -kml KML -data SIMPLIFIED_DATA_DIRECTORY``

Final Mode. This is just a wrapper around Run Mode (``-r``), which it runs for each described pair. How the numbering of pairs is done is described below under the Running Program paragraph. ``-nitfdir NITF_DIRECTORY``, ``-kml KML`` and ``-data SIMPLIFIED_DATA_DIRECTORY`` function the same way as in Run Mode. The only difference is with OUTPUT_PREFIX, since in addition 4-digit zero-padded ``ID_X`` is added to it. For example ID 127 will produce point cloud in results/OUTPUT_PREFIX0127.txt file.


Example usage:

``./mvs -f 1 tmp 1 8 12 14 16 29 30 32 34 37 50 61 88 89 108 116 "+--subpixel-mode 2 --corr-kernel 11 11 --filter-mode 2 --rm-threshold 1 --prefilter-kernel-width 1.0" -nitfdir NITF/ -kml testing/MasterProvisional1/MasterProvisional1.kml -data testing/MasterProvisional1/``

``./mvs -m 1 results/tmp????.txt provisional1.txt -kml testing/MasterProvisional1/MasterProvisional1.kml -lima 0.8 -limb 6.0 -limc 1.0 -check 10 -fill -exp -maug 0.10 -opt``

``./mvs -v provisional1.txt -vscale 2.0``

---

# Program Description

## Summary

Solution is built around NASA Ames Stereo Pipeline (https://ti.arc.nasa.gov/tech/asr/intelligent-robotics/ngt/stereo/), which produces very good point clouds from satellite imagery, but unfortunately its core functionality is limited to stereogrammetry. I use Stereo Pipeline in order to obtain point clouds from promising pairs of images. After that I use a rather simple algorithm to merge all of the point clouds which directly produces final result.


## High-Level Overview

1) Run Stereo Pipeline "stereo" executable file in order to obtain point clouds for each promising pair of NITF files. Parameters are passed by creating new stereo.default file before each run. Solution uses provided NITF files, by explicitly cropping (using "left-image-crop-win" and "right-image-crop-win" parameters) our region of interest. 

2) Convert generated point clouds to LAS files and then to text files using las2txt utility from libLAS package. Convert text files, which are in lat/long format, to UTM and crop them (with some spare margin) to required KML file.

3) Merge all point clouds using algorithm described below.

4) Optionally, post-process result in order to exploit scoring function (described below as well).


## Merging Algorithm

The first step is to compute the offsets of each image. It turns out that different pairs of images can produce results of very similar structure but up to around 10 meters apart from each other. In order to find offsets of all point clouds we use a naive greedy algorithm. We assume that the first offset is (0.0, 0.0). Then, iteratively we add next point cloud at offset that minimizes the average error between existing set of point clouds and new one. After adding all of the point clouds, we further improve the offsets by repeating the process one more time. During the second iteration we only adjust the offsets of existing point clouds (by removing them and adding them again).

One small detail that is worth mentioning is that when doing any calculation on point clouds, we convert them to 2D grid, where each tile is 30cm x 30cm. While this may lose small amount of information, effect is pretty much negligible compared to the amount of noise obtained from merging over a dozen of different point clouds. Since those 2D grid-based point clouds can be viewed as grayscale 2D images (with some pixels missing), sometimes I'm refering to those point clouds simply as images and particular points in them as pixels.

The last step is to compute the final merged point cloud. Again, my final result is simply a 2D image (with possibly some pixels missing). Each pixel (height of a tile) is computed as a median of pixels. In my final submission, I slightly deviated from this and when there's more than 5 pixels available, instead of using median, I use mean of the middle 4-5 pixels. It's hard to tell if that's a better idea. The premise behind it is that mean (after dropping potential outliers) might be better height estimation for places where it's relatively easy for photogrammetry to generate results.


## Fine-Tuning Parameters

Most of the paramaters that I use in Stereo Pipeline are default, except for these:
--subpixel-mode 2
--corr-kernel 11 11
--filter-mode 2
--rm-threshold 1
--prefilter-kernel-width 1.0

Changing subpixel-mode was obvious choice. The reason for changing other parameters was that provided satellite images were quite sharp and decent quality. At the same time, they contained a lot of small-scaled features and highly-variable terrain, which most probably is not the standard usage for Stereo Pipeline. Based on this, I assumed I should reduce the amount of blur in preprocessing. At the same time, it should be advantageous to filter out way more points than usual, since we would like to limit our point clouds to contain only high-quality results.


## Exploiting the Scoring Function

Scoring function consisted of two separate components: completeness and accuracy. Both of them could be exploited (i.e. we were able to artificially improve our results without improving our solution) by careful postprocessing of returned point cloud.

Accuracy (RMSE, Root Mean Average Error) was computed over pixels where our returned point cloud gave any value. Which means that it was possible to artificially reduce RMSE (improve accuracy), by not returning values for low-confidence places. On the other hand by doing so, we reduce completeness, so the goal was to find the balance between those two.

The easiest was to accomplish this, was to not return point clouds for areas around the edges of buildings. Or to be more exact, to check if the difference between maximum and minimum expected height in close vicinity was above some fixed threshold. I believe that in total, I was able to reduce the RMSE by around 50%.


## Running Program

First step is to download and install all of the required files/libraries. Script setup.sh contains all prerequisites. Apart from installing some libraries through "apt-get" it also downloads StereoPipeline and puts it in home directory. It doesn't have to stay there, but remember that "StereoPipeline-2.5.3-2016-08-23-x86_64-Linux/bin/" should be added to PATH in order for solution to run. It's also worth mentioning that the script downloads the data Topcoder/S3 server. If you don't need, you can just comment those lines out.

Compiling is straightforward, simply run compile.sh script.

Running the solution is a little bit more complicated. First of all, you can look at the run.sh script in order to see how the final solution was obtained. Secondly, all of the tweakable parameters are passed through command line, which allowed me a bit faster iteration process when experimenting with parameters. In general, producing new solution is split in two phases: obtaining pair-wise point clouds and merging. 

You have several ways of obtaining point clouds. The most intuitive is to run program separately for each pair (it's called "runMode" in my solution). The syntax is as follows: "./x -nitfdir PATH_TO_NITF_DIRECTORY -kml PATH_TO_KML_FILE -data PATH_TO_SIMPLIFIED_DATA_DIRECTORY -r 1 PREFIX NITF_FILE_A NITF_FILE_B [Additional parameters passed to Stereo Pipeline concatenated into a single string with preceding plus (+) sign]. ALL_CAPS means parameter. This call will use a pair of files (PATH_TO_NITF_DIRECTORY/NITF_FILE_A and PATH_TO_NITF_DIRECTORY/NITF_FILE_B), run Stereo Pipeline on them and, finally, store the point cloud in ./results/PREFIX.txt file.

There is a way to produce several point clouds, but it is a little bit more convoluted. You can run the program in "finalMode" (the same way as in run.sh script) and instead of writing the names of the files you have to provide the list of ids for pairs. The order is as follows: first sort NITF files by dates, the first 49 ids are consecutive NITF files (#0 is pair of the two oldest images and #48 is the pair of the two newest images, ids #49-#96 (48 pairs) are pairs of images with one image in-between. Similarly, #97-#143 (47 pairs) are pairs of images with two images in-between. Essentially, this is just a simple wrapper around the "runMode" in order to make testing a little bit quicker.

Fortunately, merging is pretty straightforward. After obtaining all of the point clouds (they will be placed in ./results directory) run "./x -kml PATH_TO_KML_FILE -m 1 results/tmp????.txt PATH_TO_OUTPUT_FILE [optional tweakable parameters]" and this will produce PATH_TO_OUTPUT_FILE with merged point cloud. Unfortunately, all of the parameters are not easily explained with words, so it's easiest to just read the source code in order to understand what exactly they are doing. Note that "results/tmp????.txt" is just a wildcard that matches all of the point clouds file names. You can change it to a different wildcard or just specify the file names separately.


## Additional Information

- It's quite frequently mentioned in Stereo Pipeline documentation, that RPC camera model is completely inferior when compared to camera model that is normally attached to Digital Globe satellites which is described in their XML files. In one place it's mentioned that RPC model quite often produces errors up to several meters. I never had any chance to compare those two models, but I wouldn't be surprised if the results could have been drastically improved by including DigtalGlobe's XML model. That might have explain why the only pairs of images that were remotely useful were those that were rather close in time (usually only up to one month of difference). If that's the case, it's probably safe to assume that RPC model rather drastically affects the ortorectification phase. There's also a possibility that the sizes of computed offsets during merging algorithm would have been much smaller, although this shouldn't directly affect the quality of the final result.

- Due to limited amount of provided data it wasn't practical to use any machine learning (ML) methods. While it's nearly impossible to obtain even remotely comparable results on ML alone. ML could have been used to augment the results computed from photogrammetry algorithms. One of the simplest ideas would be to use Convolutional Neural Networks (CNNs) to postprocess the merged heightmap, which hopefully would remove some of the obvious anomalies (holes in roads from cars), make builds more rectangular. It could also be used to detect building boundaries. Finally, it could also be used in order to compute confidence of the result. Biggest downside of using ML methods is that models might produce bad results when encountering new types of architecture that wasn't present in training data.

- I never tried it, but it's possible that including rotations during offset computation in merging algorithm might yield some improvement. Although, I'm quite skeptical of this since most probably I would notice something by looking at the visualizations.

- As was mentioned above, my choice of NITF pairs is hardcoded and thus this solution at this current state, cannot be used for other areas. It's possible to estimate the quality of each separate point cloud by calculating average error between considered point cloud and merged point cloud. Ideally we would also want to take confidence into consideration as well when computing our evaluation function. In order to quickly generate our "promising pairs" we can perform selection on downscaled images. 

- In order to reduce running time of whole solution, it would be best to run several instances of stereo simultaneously. Even though some parts of Stereo Pipeline are multithreaded, it runs much faster with several tasks at the same time.
