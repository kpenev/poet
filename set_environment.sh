for d in poet_src/build/libs/*/shared/release; do 
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(greadlink -e $d)
done
