// Implementation of the two query methods: LSH and brute force

#include "query.h"
#include <mpi.h>
#include <vector>
#include <utility>


void queryThreadFunct(
                     unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
                    //  unordered_map<string, GEOSGeometry*>  centeredGeosWithID,
                    unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap,
                       vector<HashMap<std::string>*> hashMaps,
                     int numNeighbors, 
                     int hashLength, 
                    // vector<unsigned int>& seeds, 
                    vector<vector<unsigned int>>& seeds, 
                     int rank,
                    // std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>>& LSHoutput,
                     ResultType& LSHoutput,
                     BoundingBox& bbox,
                     vector<vector<GEOSGeometry*>>& global_MBR_grid)
{

    // Create a new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    
    double hashingStart, hashingEnd, hashingTotal = 0.0, timeTaken = 0.0;
    size_t totalPolygons = 0;
    size_t totalHashes = 0;
    double queryingStart, queryingEnd, queryingTotal = 0.0;
    // vector<pair<double, string>> nn;  // Nearest neighbors storing only unique id's for preventing memory leak
    //  set<string> visited; // Set to keep track of visited IDs
    unordered_set<string> globalVisited;

    unordered_map<string, pair<GEOSGeometry*, vector<double>>> QuerygeometrySketchMap;

 
    // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
    for (const auto& pair : centeredQueryGeosWithID) {
        const string& id = pair.first; // Unique ID of the geometry
        GEOSGeometry* geometry = pair.second; // Geometry pointer

        vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
        if (sketchResult.empty()) {
            cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
        } else {
            QuerygeometrySketchMap[id] = make_pair(geometry, sketchResult);
        }
        //  std::cout<< " - query Sketch Size: " << sketchResult.size() << std::endl;
    
    }


    // Define Grid Size
    int gridSizeX = global_MBR_grid.size();
    int gridSizeY = global_MBR_grid[0].size();
 
    int totalLSHRefinements = 0;   // Total refined across all queries in this process
    int totalPossible = 0;         // Total brute-force comparisons (query × dataset size)

    for (const auto& item : QuerygeometrySketchMap)
    {
        const string& queryID = item.first;
        GEOSGeometry* queryGeo = item.second.first;
        const vector<double>& sketch = item.second.second;
        //vector<double> querySketch = sketch(ctx, &global_MBR_grid, queryGeo); 

        // Validate Query Geometry
        if (!queryGeo) {
            cerr << "Skipping Query ID " << queryID << " - Geometry is NULL." << endl;
            continue;
        }

        if (GEOSGeomTypeId_r(ctx, queryGeo) != GEOS_POLYGON) {
            cerr << "Skipping Query ID " << queryID << " - Not a valid Polygon." << endl;
            continue;
        }

        // Extract exterior ring before hashing
        const GEOSGeometry* exteriorRing = GEOSGetExteriorRing_r(ctx, queryGeo);
        if (!exteriorRing) {
            cerr << "Error: Failed to retrieve exterior ring for Query ID " << queryID << endl;
            continue;
        }

        // Retrieve coordinate sequence
        const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
        if (!seq) {
            cerr << "Error: Failed to retrieve coordinate sequence for Query ID " << queryID << endl;
            continue;
        }

        // Get number of points
        unsigned int numPoints;
        GEOSCoordSeq_getSize_r(ctx, seq, &numPoints);
        if (numPoints < 3) {
            cerr << "Error: Query ID " << queryID << " has less than 3 points!" << endl;
            continue;
        }

        // Convert query geometry to a vector of points
        vector<Point> queryPolygon;
        for (size_t i = 0; i < numPoints; ++i) {
            double x, y;
            GEOSCoordSeq_getX_r(ctx, seq, i, &x);
            GEOSCoordSeq_getY_r(ctx, seq, i, &y);
            queryPolygon.push_back({x, y});
        }
        
     
     
        // Precompute occupied grid cells
        unordered_set<int> occupiedCells = precomputeOccupiedCells(sketch, gridSizeX, gridSizeY);

     
        // Declare nn properly as a vector of (distance, neighborID) 
        std::vector<std::pair<double, std::string>> nn;  
        //set<string> visited; // Track visited IDs
        unordered_set<string> visitedForQuery;  //  Prevents duplicates across hashmaps (but resets for each query)


        // Create a vector to store hashing times for each query polygon
        unordered_map<string, double> queryHashingTimes;

        int refinementCount = 0;
        int possible = geometrySketchMap.size();  // brute-force baseline
        //unordered_set<string> visitedForQuery;

        for (unsigned int h = 0; h < hashMaps.size(); h++) {
            hashingStart = MPI_Wtime();

            vector<int> hash = hash2d(queryPolygon, hashLength, bbox, global_MBR_grid, occupiedCells, seeds[h], ctx);

            hashingEnd = MPI_Wtime();
            hashingTotal += hashingEnd - hashingStart;

            queryingStart = MPI_Wtime();

            const list<string>* bucket = hashMaps[h]->get(hash);

            if (bucket != nullptr) {
                for (const string& neighborID : *bucket) {
                    if (visitedForQuery.find(neighborID) == visitedForQuery.end()) {
                        visitedForQuery.insert(neighborID);

                        auto neighborSketchIt = geometrySketchMap.find(neighborID);
                        if (neighborSketchIt != geometrySketchMap.end()) {
                            GEOSGeometry* neighborGeo = neighborSketchIt->second.first;
                            double distance = jaccardDistance(ctx, queryGeo, neighborGeo);
                            refinementCount++; 
                            nn.push_back(make_pair(distance, neighborID));
                        }
                    }
                }
            }

            queryingEnd = MPI_Wtime();
            queryingTotal += queryingEnd - queryingStart;
        } // End of hashMaps loop


        totalLSHRefinements += refinementCount;
        totalPossible += possible;
   
        // Sort the nearest neighbors by distance
        sort(nn.begin(), nn.end(), [](const std::pair<double, std::string>& a, const std::pair<double, std::string>& b) {
            return a.first < b.first;  // Compare distances
        });


        // Keep only the required number of nearest neighbors
        if (nn.size() > numNeighbors) {
            nn.resize(numNeighbors);
        }
        

        // Store the nearest neighbors for the current query
        LSHoutput.push_back(make_pair(queryID, nn));
    }

    int globalLSHRefinements = 0;
    int globalPossible = 0;
    double maxhashingTotal;
    double maxqueryingTotal;
    // Use MPI_Reduce to find the maximum hashing time across all processes
    MPI_Reduce(&hashingTotal, &maxhashingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Use MPI_Reduce to find the maximum hashing time across all processes
    MPI_Reduce(&queryingTotal,&maxqueryingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&totalLSHRefinements, &globalLSHRefinements, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalPossible, &globalPossible, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Max LSH Hashing Time for query geometry: " << maxhashingTotal << " seconds" << std::endl;
        std::cout << "Max Querying Time (hash map search for neighbor and refinement): " << maxqueryingTotal << " seconds" << std::endl;

        std::cout << "\n=== Global LSH Pruning Stats ===" << std::endl;
        std::cout << "Total refined (across all queries): " << globalLSHRefinements << std::endl;
        std::cout << "Total brute-force comparisons: " << globalPossible << std::endl;
        std::cout << "Overall pruning ratio: " 
              << (1.0 - (double)globalLSHRefinements / globalPossible) * 100.0 << "%" << std::endl;
    }

    // Clean up the GEOS context
    GEOS_finish_r(ctx);
     
}



//copy geometries + id's
vector<GEOSGeometry *> *copyGeometriesFromMap(const unordered_map<string, GEOSGeometry *> *original, GEOSContextHandle_t ctx) 
{
    vector<GEOSGeometry *> *copy = new vector<GEOSGeometry *>();
    for (const auto &pair : *original) {
        // Make a copy of the geometry
        GEOSGeometry *copyGeometry = GEOSGeom_clone_r(ctx, pair.second);
        copy->push_back(copyGeometry);
    }
    return copy;
}



// Function to print the contents of each HashMap in the vector
void printAllHashMaps(const std::vector<HashMap<std::string>*>& maps) {
    for (size_t i = 0; i < maps.size(); ++i) {
        std::cout << "Hash Map " << i << ":" << std::endl;
        // Call the print method on the current HashMap instance
        maps[i]->print(true); 
    }
}

//generate lsh hashes of geos geometries for further querying
void LSHHashFunction(
    const unordered_map<string, GEOSGeometry*>& geos, // Changed to unordered_map
    int rank,
    const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
    unsigned int numNeighbor,
    int gridSize,
    int hashLength ,
   // std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>>& LSHoutput,
    ResultType& LSHoutput,
    GEOSContextHandle_t ctx) 
{

 // Create a new GEOS context handler
   
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    GEOSWKTWriter* wktWriter = GEOSWKTWriter_create_r(ctx);
    
    //copy original geometries as when original geometries was used, it was resulting in segmentation fault
    unordered_map<string, GEOSGeometry*> geosCopy = copyGeometriesMap(geos, ctx);
  
    // Vector to hold centered geometries
    unordered_map<string, GEOSGeometry*> centeredGeosWithID;

    // Iterate over the unordered_map, center each geometry, and store it with its unique id
    for (const auto& pair : geos) {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            // Store the centered geometry with its identifier in the map
            centeredGeosWithID[pair.first] = centeredGeo;
        } else {
            std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
        }
    }

    vector<GEOSGeometry*> geometries;
    for (const auto& pair : geosCopy) {
        GEOSGeometry* originalGeometry = pair.second;
        // Call the centerGeometry function on the original geometry
        GEOSGeometry* centeredGeometry = centerGeometry(ctx, originalGeometry);
        if (centeredGeometry != nullptr) {
            // Successfully centered the geometry, add it to the vector
            geometries.push_back(centeredGeometry);
            GEOSGeom_destroy_r(ctx, originalGeometry);
        } else {    
            // Centering failed, handle the error accordingly
            cerr << "Warning: Failed to center geometry for ID " << pair.first << endl;
        }
    }

    double hashingStart, hashingEnd, hashingTotal = 0.0;
    double insertStart, insertEnd, insertTotal = 0.0;

    // Start timing LSH
    double startTimeLSH = MPI_Wtime();
    
    // Start construction timing  for polygon MBR, grids, sketch 
    double startTimeConstruction = MPI_Wtime();

    // Start construction timing  for polygon MBR, grids, sketch 
    double startTimeMBR = MPI_Wtime();
 
    //Create Global MBR parallely by calling GlobalMbr_parallel function
    GEOSGeometry *global_MBR = GlobalMbr_parallel(&geometries, rank);

    //Creating grid over MBR in root only
    vector<vector<GEOSGeometry *>> global_MBR_grid = createGrid(ctx,
                                            global_MBR, gridSize); 
                                            
    unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap;
                                            
    // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
    for (const auto& pair : centeredGeosWithID) 
    {
        const string& id = pair.first; // Unique ID of the geometry
        GEOSGeometry* geometry = pair.second; // Geometry pointer

        vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
        if (sketchResult.empty()) 
        {
            cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
        } else 
        {
            geometrySketchMap[id] = make_pair(geometry, sketchResult);
        }
    }


    // OVERALL time for construction 
    double endTimeMBR = MPI_Wtime();
    double totalTimeMBR = endTimeMBR - startTimeMBR;
    
    double maxTotalTimeMBR;

    // Use MPI_Reduce to find the maximum hashing time across all processes
    MPI_Reduce(&totalTimeMBR, &maxTotalTimeMBR, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    // The master process (rank 0) prints the maximum hashing time
    if (rank == 0)
    {
        std::cout << "time for construction of globalMBR for input: " << totalTimeMBR << " seconds" << std::endl;
    }

    int nmaps = 2;

    std::vector<HashMap<string>*> maps(nmaps);
    for (size_t i = 0; i < nmaps; ++i)
    {
        maps[i] = new HashMap<string>(10000, 1, nullptr); 
    } // Adjust the template parameter to int

    //  Generate different seeds per hash table
    vector<unsigned int> baseSeeds = {100};  // Base seed set
    vector<vector<unsigned int>> tableSeeds(nmaps, vector<unsigned int>(hashLength));


    // Get the envelope of the global MBR
    GEOSGeometry* mbrEnvelope = GEOSEnvelope_r(ctx, global_MBR);
    if (!mbrEnvelope) {
        cerr << "Error: Failed to compute envelope for the global MBR." << endl;
        return;
    }


    BoundingBox globalBoundingBox = extractBoundingBox(ctx, mbrEnvelope);
  
    
    for (const auto& pair : geometrySketchMap) {
        const string& id = pair.first;
        GEOSGeometry* centeredGeometry = pair.second.first; // Get centered geometry
        const vector<double>& sketch = pair.second.second;  // Get associated sketch vector
    
        // Define Grid Size
        int gridSizeX = global_MBR_grid.size();
        int gridSizeY = global_MBR_grid[0].size();
    
        // Precompute occupied grid cells
        unordered_set<int> occupiedCells = precomputeOccupiedCells(sketch, gridSizeX, gridSizeY);

        // Validate before calling hash2d
        if (!centeredGeometry) {
            cerr << "Skipping ID " << id << " - Geometry is NULL." << endl;
            continue;
        }
    
        if (GEOSGeomTypeId_r(ctx, centeredGeometry) != GEOS_POLYGON) {
            cerr << "Skipping ID " << id << " - Not a valid Polygon." << endl;
            continue;
        }

        // Extract exterior ring
        const GEOSGeometry* exteriorRing = GEOSGetExteriorRing_r(ctx, centeredGeometry);
        if (!exteriorRing) {
            cerr << "Error: Failed to retrieve exterior ring for polygon ID " << id << endl;
            continue;
        }
    
        // Retrieve coordinate sequence
        const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
        if (!seq) {
            cerr << "Error: Failed to retrieve coordinate sequence for polygon ID " << id << endl;
            continue;
        }

        // Get number of points
        unsigned int numPoints;
        GEOSCoordSeq_getSize_r(ctx, seq, &numPoints);
    
        if (numPoints < 3) {
            cerr << "Error: Polygon ID " << id << " has less than 3 points!" << endl;
            continue;
        }

        // Convert to vector of points
        vector<Point> polygon;
        for (size_t i = 0; i < numPoints; ++i) {
            double x, y;
            GEOSCoordSeq_getX_r(ctx, seq, i, &x);
            GEOSCoordSeq_getY_r(ctx, seq, i, &y);
            polygon.push_back({x, y});
        }

        try {
            hashingStart = MPI_Wtime();
 
            for (size_t h = 0; h < nmaps; ++h) {
        
                for (size_t i = 0; i < hashLength; ++i) {
                    tableSeeds[h][i] = baseSeeds[i % baseSeeds.size()] + h * 123 + i * 17;  // Different per table & per seed
                }
    
                vector<int> hash = hash2d(polygon, hashLength, globalBoundingBox, global_MBR_grid, occupiedCells, tableSeeds[h], ctx);

                maps[h]->insert(hash, id);
            
                // Compute area of polygon
                // Compute area of polygon
                double polygonArea = 0.0;
                if (!GEOSArea_r(ctx, centeredGeometry, &polygonArea)) {
                    std::cerr << "Warning: Failed to compute area for polygon ID " << id << std::endl;
                    continue;
                }

                // Compute area of global MBR (reuse bounding box)
                double globalMBRArea = (globalBoundingBox.xMax - globalBoundingBox.xMin) *
                                       (globalBoundingBox.yMax - globalBoundingBox.yMin);

                if (globalMBRArea == 0) {
                    std::cerr << "Error: Global MBR has zero area." << std::endl;
                    continue;
                }

                double sparsity = polygonArea / globalMBRArea;

                // Estimate total bits needed to store this hash vector
                int totalBits = 0;
                for (int val : hash) {
                    totalBits += (val > 0) ? static_cast<int>(std::ceil(std::log2(val + 1))) : 1;
                }

                // Only rank 0 writes to file
                if (rank == 0) {
                    std::ofstream bitOut("/home/asbmr/LSH_SIM_CEM_2d_multiplemaps/similarity_search_project/geometric-ANN-main/data/polygon_sparsity_test_urban_hashbits.csv", std::ios::app);
                    if (bitOut.is_open()) {
                        bitOut << id << "," << sparsity << "," << totalBits << "\n";
                        bitOut.close();
                    } else {
                        std::cerr << "Rank 0 could not open CSV file for writing.\n";
                    }
                }

            }

            hashingEnd = MPI_Wtime();
            hashingTotal += hashingEnd - hashingStart;


        } 
        catch (const std::exception& e) {
            std::cerr << "Error hashing sketch for ID " << id << ": " << e.what() << std::endl;
        }
    }

    //printAllHashMaps(maps);   

    //for hshing and hashmap insertion 
    double maxhashingTotal;
    double maxinsertTotal;

    // Use MPI_Reduce to find the maximum hashing time across all processes
    MPI_Reduce(&hashingTotal, &maxhashingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&insertTotal, &maxinsertTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "LSH Hashing Time: " << maxhashingTotal << " seconds" << std::endl;
        std::cout << "hash map insertion Time: " << maxinsertTotal << " seconds" << std::endl;
    }
 
    queryThreadFunct(centeredQueryGeosWithID, geometrySketchMap, maps, numNeighbor, hashLength, tableSeeds, rank, LSHoutput, globalBoundingBox, global_MBR_grid);

    if (LSHoutput.empty()) {
        std::cout << "LSHoutput is empty." << std::endl;
    }


    // Cleanup
    for (auto& pair : centeredGeosWithID) {
        GEOSGeom_destroy_r(ctx, pair.second);
    }
    //}
    GEOSGeom_destroy_r(ctx, global_MBR); 
}


