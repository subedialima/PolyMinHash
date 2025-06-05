# PolyMinHash
PolyMinHash: Efficient Area-Based MinHashing of Polygons for Approximate Nearest Neighbor Search

This repository contains the source code for PolyMinHash.

This research hopes to use **rejection sampling** to find the nearest neighbors of polygons using Jaccard Distance as our metric. The project is primarily written in C++ using the reentrant C API for the GEOS library for the hashing algorithms. Data visualization is done with Python using the Shapely along with matplotlib.

The project currently has multiple directories.

data Where the input data is stored and the output data is written. Due to the large nature of the datasets, data is not available in this repository but can be found at the UCR Star website.
src The actual source code for the nearest neighbor search code, written in C++20.
Within the code itself, each function, class, type, etc in a header file has a description of its respective parameters, return values, members, etc, which I hope will make understanding the code easier.


