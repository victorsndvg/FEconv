(copiado de almacen/code/basicmod)
dependencies on basicmod modules:

    compiler, os -> report -> convers -> files

dependencies on basicmod/alloc modules:

    ... -> report -> alloc \
                            > system
                     files / 
    
In summary:

    compiler, os -> report -> convers, alloc 
                                 |          \
                                 v           > system
                               files -------/
    
