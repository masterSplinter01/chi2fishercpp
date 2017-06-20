#!/bin/bash

read -e -p "Enter dir path, use tab for completion: " srcdir
read -e -p "Where do u want to place the package?: " resdir
cd $resdir

R CMD check $srcdir

if [ $? -eq 0 ]; then
    echo BUILDING PACKAGE...
    R CMD build $srcdir
else
	echo FAIL
fi


