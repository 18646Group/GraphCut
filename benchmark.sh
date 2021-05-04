BINARY_PATH=$1

for file in brick.gif green.gif greenLeaf.jpg
do
    echo "processing file $file"
    $BINARY_PATH images/originals/$file images/outputs/$file 512x512 | grep "Run time"
done
