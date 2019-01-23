#!/bin/sh

sfroot=$PWD

cd examples
dirs=()
ndir=0
for FILE in `ls -l`
do
    if test -d $FILE; then
      dirs+=( $FILE )
      ndir=$[ndir+1]
    fi
done

echo "Enter the index of the example to run (1-$ndir):"

for ((i=0; i<$ndir; i++)); do      
    echo "$[i+1]) ${dirs[i]}"
done

while true; do
  read ie
  re='^[0-9]+$'
  if ! [[ $ie =~ $re ]]; then
    echo "Please input a valid number"
  elif [ $ie -gt 0 ] && [ $ie -le $ndir ]; then
    break
  else
    echo "Please input a number between 1 and $ndir"
  fi
done

dir=${dirs[ie-1]}
echo "Running example/"$dir

cd $dir

rm -rf output*
rm -rf scratch

export PYTHONPATH=$sfroot
python $sfroot"/scripts/sfrun" --workdir=$PWD