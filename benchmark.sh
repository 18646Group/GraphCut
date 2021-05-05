#Usage: bash benchmark.sh ./graph_cut_omp 256x256

BINARY_PATH=$1
OUTPUT_SIZE=$2

for file in brick.gif green.gif greenLeaf.jpg MoreRocks.gif akeyboard_small.gif larva.gif t.yello.01.gif
do
    echo "processing file $file into $OUTPUT_SIZE"
    $BINARY_PATH images/originals/$file images/outputs/$file $OUTPUT_SIZE | grep "Run time"
done
