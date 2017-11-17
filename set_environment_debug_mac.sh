for d in build/libs/*/shared/debug; do 
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(greadlink -e $d)
done
