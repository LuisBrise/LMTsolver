#!/bin/bash

# Function to check if input is an integer
is_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

# Prompt user for an integer input
read -p "Enter the number of folders to create in each directory: " n

# Validate the input
if is_integer "$n"; then
    # Loop through all folders in the current directory
    for dir in */; do
        if [ -d "$dir" ]; then
            echo "Creating $n folders in $dir"
            for i in $(seq 1 "$n"); do
                mkdir "$dir/folder_$i"
            done
        fi
    done
    echo "Folders have been created in all subdirectories."
else
    echo "Error: '$n' is not a valid integer."
    exit 1
fi
