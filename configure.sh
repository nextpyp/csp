#!/bin/bash

if ! [ -n "$1" ]
then
    echo "Usage: ./configure.sh. /PATH/TO/FREALIGN/EXECUTABLE/frealign_v8.exe"
    exit
else
    frealign_executable=$1
fi

# FOLDERS='bin data frealign/log frealign/maps frealign/scratch frealign/swarm'
FOLDERS='bin'
for k in $FOLDERS
do 
    echo -ne "Creating folder $k: "
    mkdir -p $k
    if [ "$?" -eq "0" ]; then
        echo "[DONE]"
    else
        echo "[FAIL]"
        exit
    fi
done

echo "Running CMAKE:"
cd bin
cmake ../src/
cd ..

ln -s -f $frealign_executable frealign/frealign_v8.exe

echo -e "\n=> To build CSP run ./make_csp.sh"
