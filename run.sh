#!/bin/sh
export PATH=$PATH:~/StereoPipeline-2.5.3-2016-08-23-x86_64-Linux/bin/

echo $PATH

KML="$1"
NITF="$2"
DATA="$3"
OUTPUT="$4"

./x -nitfdir $NITF -kml $KML -data $DATA -output $OUTPUT -f 1 tmp 1 8 12 14 16 29 30 32 34 37 50 61 88 89 108 116 "+--subpixel-mode 2 --corr-kernel 11 11 --filter-mode 2 --rm-threshold 1 --prefilter-kernel-width 1.0"
./x -kml $KML -m 1 results/tmp????.txt $OUTPUT -lima 0.8 -limb 6.0 -limc 1.0 -check 10 -fill -exp -maug 0.10 -opt
rm results/*