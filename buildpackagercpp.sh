#!/bin/bash

read -e -p "Enter filename, use tab for completion: " file

R CMD check $file

if [ $? -eq 0 ]; then
    echo BUILDING PACKAGE...
    R CMD build $file
else
	echo FAIL
fi


