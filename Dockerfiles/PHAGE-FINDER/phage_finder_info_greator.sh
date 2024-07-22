#!/bin/bash
# This script performs file manipulation to create the phagefinder input file = phage_finder_info.txt
# Usage . phage_finder_info_generator.sh <annotation gb file>

# User Input and Arguments
show_help() {
    echo "Usage: $0 [-o operation] <file>"
    echo "Options:"
    echo "  -o operation    Specify the operation to perform (e.g., sort, unique)."
    echo "  -v              Enable verbose mode."
    echo "  -h              Show this help message"
}

# User Feedback
verbose = false

while getopts ":o:hv" opt; do 
    case ${opt} in 
        o )
            operation=$OPTARG
            ;;
        v )
            verbose=true
            ;;
        h )
            show_help
            exit )
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            show_helpexit1
            ;;
        : )
            echo "Invalid option: -$OPTARG requires and argument" 1>&2
            show_helpexit 1
            ;;
    esac
done

done shift $((OPTIND -1))
file=$1

# Input Validation & Error Handling
if [ -z "$operation" ]; then
    echo "Operation is required."
    show_help
    exit 1
fi

if [ -z "$file" ]; then
    echo "File is required."
    show_help
    exit 1
fi

if [ ! -f "$file" ]; then
    echo "File not found: $file"
    exit 1
fi

if [ ! -r "$file" ]; then
    echo "File is not readable: $file"
    exit 1
fi

if [ "$verbose" = true ]; then
    echo "Operation: $operation"
    echo "File: $file"
fi

#Performing Operations
case $operation in
    sort)
        sort "$file"
        ;;
    unique)
        sort "$file" | uniq
        ;;
    *)
        echo "Unknown operation: $operation"
        show_help
        exit 1
        ;;
esac
