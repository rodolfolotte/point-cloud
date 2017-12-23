#! /bin/sh

# Example
#./preprocessing PATH/$1.$2 PATH/$1-simplified.$2 -simplify wlop -nb_wlop_simplify 0.5 -retain_percentage 10.0
#./preprocessing PATH/$1-simplified.$2 PATH/$1-normal.$2 -normal plane -nb_neighbors_pca 18 
#./preprocessing PATH/$1-normal.$2 PATH/$1-oriented.$2 -orient -nb_neighbors_mst 18 
#./preprocessing PATH/$1-oriented.$2 PATH/$1-structured.$2 -structure -epsilon 0.015 -paral true -ortho true -coplan false -zsymmetry true -degree 10

./preprocessing $1.$2 $1-simplified.$2 -simplify wlop -nb_wlop_simplify 0.5 -retain_percentage 10.0
./preprocessing $1-simplified.$2 $1-normal.$2 -normal plane -nb_neighbors_pca 18 
./preprocessing $1-normal.$2 $1-oriented.$2 -orient -nb_neighbors_mst 18 
./preprocessing $1-oriented.$2 $1-structured.$2 -structure -epsilon 0.015 -paral true -ortho true -coplan false -zsymmetry true -degree 10