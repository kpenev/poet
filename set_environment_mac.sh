for d in build/libs/*/shared/release; do 
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(greadlink -e $d)
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(greadlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(greadlink -e $d)
done
