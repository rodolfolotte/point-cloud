# Command-line point-cloud preprocessing
This is a collection CGAL algorithms, which I have selected the main procedures for preprocessing point-clouds. The code is an adaptation using different examples and demo from CGAL initiative, wherein the adaptations were thinking to be an easy to use command-line for basic point-cloud processing. Such as:

Possible file extention: ply, obj, off, pwd, and xyz.
Details about the method and implementation: https://doc.cgal.org/latest/Point_set_processing_3/index.html.

## SIMPLIFICATION
Possible simplification methods: wlop, random, hierarchy, and var_max
Possible flags: -nb_wlop_simplify, -retain_percentage

Usage: `./preprocessing INPUT-POINT-CLOUD OUTPUT-POINT-CLOUD -simplify METHOD -nb_wlop_simplify FLOAT_VALUE -retain_percentage FLOAT_VALUE`
Example: `./preprocessing samples/fandisk.off samples/fandisk-simplified.off -simplify wlop -nb_wlop_simplify 0.5 -retain_percentage 10.0`

## SMOOTH
**NOT-IMPLEMENTED**

## OUTLIER REMOVAL
**NOT-IMPLEMENTED**

## NORMAL ORIENT ESTIMATION
Possible estimate normal methods: plane, quadric, and vcm (default=plane)
Possible flags: 
	-nb_neighbors_pca: Number of neighbors to compute tangent plane (default=18)"
    -nb_neighbors_jet_fitting: Number of neighbors to compute quadric (default=18)
    -offset_radius_vcm: Offset radius to compute VCM (default=0.1)
    -convolve_radius_vcm: Convolve radius to compute VCM (default=0.1)

Usage: `./preprocessing INPUT-POINT-CLOUD OUTPUT-POINT-CLOUD -normal METHOD FLAG XXX_VALUE`
Example: `./preprocessing samples/fandisk.off samples/fandisk-normals.off -normal -nb_neighbors_pca 18`

## ORIENT NORMALS
Possible flags: 
	-nb_neighbors_mst: Number of neighbors to compute the MST (default=18)

Usage: `./preprocessing INPUT-NORMAL-POINT-CLOUD OUTPUT-ORIENTED-POINT-CLOUD -orient FLAG INT_VALUE`
Example: `./preprocessing samples/fandisk-normals.off samples/fandisk-oriented.off -orient -nb_neighbors_mst 12`

## STRUCTURE
Possible flags: 
	-epsilon: (default=0.015)
    -paral: Regularize parallelism? (default=true)
    -ortho: Regularize orthogonality? (default=true)
    -coplan: Regularize coplanarity? (default=false)
    -zsymmetry: Regularize z-symmetry? (default=true)
    -degree: Degree of tolerance (default=10)"

Usage: `./preprocessing INPUT-ORIENTED-POINT-CLOUD OUTPUT-STRUCTURED-POINT-CLOUD -structure FLAGS VALUES`
Example: `./preprocessing samples/fandisk-oriented.off samples/fandisk-structured.off -structure -epsilon 0.015 -paral true -ortho true -coplan false -zsymmetry true -degree 10`
    
## MANY DATASETS
To be used in a procedural way, with different sets of data sequentially, a shell-script (run.sh) is also provided.

Usage: `# Usage: ./run.sh POINT_CLOUD_FILE_WITHOUT_EXTENSION FILE_EXTENSION_WITHOUT_DOT`

# Refactoring
The code was not implemented for high performance or to meet any development metrics. Therefore, it is expected a refactoring in which will be considered the proper implementation of classes and funtionality.