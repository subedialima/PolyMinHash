# PolyMinHash

**PolyMinHash: Efficient Area-Based MinHashing of Polygons for Approximate Nearest Neighbor Search**

This repository contains the source code for the PolyMinHash system, which performs approximate nearest neighbor (ANN) search for polygon datasets using an area-based rejection sampling approach. The method uses MinHash-style hashing combined with polygon geometry to enable fast, scalable similarity search based on Jaccard distance.

---

## 📌 Overview

PolyMinHash uses rejection sampling to generate hash values based on polygon area sparsity. Instead of traditional grid encodings or text-based representations, our method encodes each polygon as a compact hash vector, enabling efficient filtering and ANN queries.

* **Hashing Method**: Rejection sampling based on 2D point-in-polygon tests.
* **Similarity Metric**: Jaccard distance based on geometric intersection area.
* **Pruning**: Achieves up to 98% pruning while maintaining recall of 69% to 97%.
---

## 📁 Repository Structure

```
PolyMinHash/
├── src/                 # C++20 source code
│   ├── main.cpp         # Entry point
│   ├── mpi_gis.cpp      # MPI parallelism
│   ├── geoutil.cpp      # Geometry utilities
│   └── query.cpp        # Minhashing
│   └── brute_force.cpp  # Brute-Force implementation(exact method)
│   └── util.h.cpp       # Basic utilities (other)
│   └── parse_geodata.cpp# read wkt files
├── data/                # Input WKT files and output CSVs (not included)
├── README.md            # This file

---

## ⚙️ Dependencies

* GEOS (C API) — tested with GEOS 3.12+
* MPI (e.g., OpenMPI)
* C++20 compiler (e.g., g++ or clang++)
* Python (optional for visualization)

  * `shapely`, `matplotlib`

---

## 🔧 Building

```bash
make
```

---

## 🚀 Running

Example command:

```bash
mpirun -n 100 ./spjoin data/input_polygon.wkt data/query_polygon.wkt 10

```

---


## 📦 Datasets

Due to size constraints, input datasets are not hosted in this repository. Publicly available WKT-format polygon datasets (e.g., from the UCR STAR repository) should be placed in the `data/` directory.

---


## 📬 Contact

For questions or collaborations, contact: [asbmr@mst.edu](mailto:asbmr@mst.edu)  [subedialima0@gmail.com](mailto:subedialima0@gmail.com)

---
