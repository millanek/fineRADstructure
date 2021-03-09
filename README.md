#  fineRADstructure
Publication:  
Malinsky, M., Trucchi, E., Lawson, D. J. & Falush, D. (2018) RADpainter and fineRADstructure: Population Inference from RADseq Data. Molecular Biology and Evolution 35, 1284–1290. doi: [https://doi.org/10.1093/molbev/msy023](https://doi.org/10.1093/molbev/msy023 ) 

## Quickstart:
### Requirements:
Assuming that you are in a Linux/UNIX like environment, you need:
1. a) On Linux/Unix a reasonably newish version of the gcc compiler (gcc 4.6.2 and onwards work fine)
1. b) Mac/Apple users need to have the Command Line Tools installed
2. GNU Scientific Development Library, (libgsl0-dev and libgsl0 in Ubuntu, gsl-devel in OpenSUSE) 

### Installation:
Briefly:
`./configure` 
`make`  

This package includes two exacutables: `RADpainter` and `fineSTRUCTURE`

### Usage:
1. Calculate the co-ancestry matrix:  `./RADpainter paint INPUT_RAD_FILE.txt`
2. Assign individuals to populations:  `./finestructure -x 100000 -y 100000 -z 1000 INPUT_RAD_FILE_chunks.out INPUT_RAD_FILE_chunks.mcmc.xml`
3. Tree building:  `./finestructure -m T -x 10000 INPUT_RAD_FILE_chunks.out INPUT_RAD_FILE_chunks.mcmc.xml INPUT_RAD_FILE_chunks.mcmcTree.xml`
 
## Further information:

Detailed instructions are available on my [personal homepage](https://www.milan-malinsky.org/fineradstructure). Please also read the publication.

`RADpainter` was written by Milan Malinsky (millanek@gmail.com) and `fineSTRUCTURE` was written by Daniel Lawson (dan.lawson@bristol.ac.uk)
