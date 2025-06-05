#include "geoutil.h"
#include <mpi.h>

vector<GEOSGeometry *> rootAggregateMBRs;
GEOSGeometry *boundingBoxGeometry = nullptr;
GEOSGeometry *globalMBR = nullptr;


//to gather MBR coordinates at root from each processes
std::vector<std::pair<double, double>> gatherCoordinatePairs(const std::vector<std::pair<double, double>>& coordinates, int rank)
{
    std::vector<std::pair<double, double>> allCoordinates;

    // Calculate the total number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate space on the root process to store gathered values
    if (rank == 0)
    {
        allCoordinates.resize(size * coordinates.size());
    }

    // Gather values from all processes to the root process
    MPI_Gather(coordinates.data(), 2 * static_cast<int>(coordinates.size()), MPI_DOUBLE,
               allCoordinates.data(), 2 * static_cast<int>(coordinates.size()), MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Return the gathered coordinates to all processes
    return allCoordinates;
}




//To create Global MBR at root and then broadcast it to every other processes
GEOSGeometry  *GlobalMbr_parallel(vector<GEOSGeometry *> *geos, int rank)
{ 
// Create new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    
    GEOSGeometry *aggregateLocalMBR = nullptr;
    double localMBRStart = MPI_Wtime();  // Start timer for local MBRs
    vector<GEOSGeometry *> localMBRs;

    // Calculate MBRs for each polygon in the specified range
   for (std::vector<GEOSGeom_t*>::size_type i = 0; i < geos->size(); ++i)
    {
        GEOSGeometry *currentPolygon = geos->at(i);
      
        GEOSGeometry *mbr = minimumBoundingRectangle(ctx, currentPolygon);
        if (!mbr)
        {
            cerr << "Error: Null MBR for polygon " << i << " in process " << rank << endl;
            continue;
        }
        localMBRs.push_back(mbr); // Store the local MBR in the vector
        
        GEOSGeom_destroy_r(ctx, currentPolygon);
       
    }
    
    double localMBREnd = MPI_Wtime();
    double localMBRTime = localMBREnd - localMBRStart;
    
    if (rank == 0) {
    std::cout << "[Rank 0] Time to compute local MBRs: in geoutil.cpp " 
              << localMBRTime << " seconds" << std::endl;
   }

    // Calculate aggregate local MBR for each process
    GEOSGeometry *LocalMBR_union = calculateAggregateLocalMBR(ctx, localMBRs);
     
    
    aggregateLocalMBR = minimumBoundingRectangle(ctx, LocalMBR_union);
    
    GEOSGeom_destroy_r(ctx, LocalMBR_union);
   
    //convert mbr coordinates to double type 
    std::vector<std::pair<double, double>> mbrCoordinates = extractCoordinates(ctx, aggregateLocalMBR);
  
    //gather all double coordinates to root
    std::vector<std::pair<double, double>> mbrCoordinates_gathered  = gatherCoordinatePairs(mbrCoordinates,  rank);
      
    //extracting only MBR coordinates
     auto boundingBox = findBoundingBox(mbrCoordinates_gathered);
     
    //broadcasting double type MBR coordinates to all process
    MPI_Bcast(&boundingBox, sizeof(std::pair<std::pair<double, double>, std::pair<double, double>>), MPI_BYTE, 0, MPI_COMM_WORLD);
 
    //creating geos polygon from double data type
    boundingBoxGeometry = createPolygonFromBoundingBox(ctx, boundingBox);
     
    globalMBR = minimumBoundingRectangle(ctx, boundingBoxGeometry);
   
    GEOSGeom_destroy_r(ctx, boundingBoxGeometry);

return globalMBR;

}
    

//to copy geometries with it's unique id 
unordered_map<string, GEOSGeometry*> copyGeometriesMap(
const unordered_map<string, GEOSGeometry*>& originalMap, GEOSContextHandle_t ctx
) 
{
    unordered_map<string, GEOSGeometry*> copiedMap;
    for (const auto& pair : originalMap) {
        const string& id = pair.first;
        GEOSGeometry* originalGeometry = pair.second;
        // Clone the geometry using GEOSGeom_clone_r
        GEOSGeometry* copiedGeometry = GEOSGeom_clone_r(ctx, originalGeometry);
        // Insert the id and copied geometry into the new map
        copiedMap.insert({id, copiedGeometry});
    }
    return copiedMap;
}


//to calculate Jaccard distance for comparison
double jaccardDistance(GEOSContextHandle_t ctx, GEOSGeometry *g1,
                       GEOSGeometry *g2)
{
    GEOSGeometry *Intersect = GEOSIntersection_r(ctx, g1, g2);
    double intersectArea;
    GEOSArea_r(ctx, Intersect, &intersectArea);
    GEOSGeom_destroy_r(ctx, Intersect);
    
    // Print intersection area
   // std::cout << "Intersection area: " << intersectArea << std::endl;

    GEOSGeometry *Union = GEOSUnion_r(ctx, g1, g2);
    double unionArea;
    GEOSArea_r(ctx, Union, &unionArea);
    GEOSGeom_destroy_r(ctx, Union);
    
    // Print union area
   // std::cout << "Union area: " << unionArea << std::endl;

    return 1.0 - (intersectArea / unionArea);
}

double SketchJaccardDistance(const std::vector<double>& sketch1, const std::vector<double>& sketch2) {
    if (sketch1.size() != sketch2.size()) {
        throw std::invalid_argument("Sketches must be of the same size.");
    }

    double minSum = 0.0; // Sum of minimum values for each cell
    double maxSum = 0.0; // Sum of maximum values for each cell

    for (size_t i = 0; i < sketch1.size(); ++i) {
        minSum += std::min(sketch1[i], sketch2[i]);
        maxSum += std::max(sketch1[i], sketch2[i]);
    }

    // Avoid division by zero in case maxSum is zero
    if (maxSum == 0) {
        return 1.0; // If both sketches are empty, consider them completely dissimilar
    }

    double jaccardSimilarity = minSum / maxSum;

    // Jaccard distance is 1 minus Jaccard similarity
    return 1.0 - jaccardSimilarity;
}

double SketchJoinDistance(const std::vector<double>& sketch1, const std::vector<double>& sketch2) {
    if (sketch1.size() != sketch2.size()) {
        throw std::invalid_argument("Sketches must be of the same size.");
    }

    double intersection = 0.0; // Sum of the products of corresponding elements
    double unionSum = 0.0;     // Sum of all elements in both sketches

    for (size_t i = 0; i < sketch1.size(); ++i) {
        intersection += sketch1[i] * sketch2[i]; // Calculate intersection
        unionSum += sketch1[i] + sketch2[i];     // Calculate union
    }

    // Union is unionSum minus the intersection
    unionSum -= intersection;
    
  // Print the intersection area
   // std::cout << "Intersection area of 2 sketches: " << intersection << std::endl;

    

    // Avoid division by zero in case unionSum is zero
    if (unionSum == 0) {
        return 1.0; // If both sketches have no elements, consider them completely dissimilar
    }

// Print the union area
   // std::cout << "Union area of sketches: " << unionSum << std::endl;

    double jaccardDistance = intersection / unionSum;

    // Return the Jaccard distance
    return 1.0 - jaccardDistance; // Jaccard distance is 1 minus Jaccard similarity
}


double fillRatio(GEOSContextHandle_t ctx, GEOSGeometry *base,
                 GEOSGeometry *overlay)
{  
 //check for validity
  int isValidb = GEOSisValid_r(ctx, base);
    if (GEOSisValid_r(ctx, base) != 1)
        aprintf(28, "Invalid Grid Cell detected!\n");
     
     int isValido = GEOSisValid_r(ctx, overlay);
    if (GEOSisValid_r(ctx, overlay) != 1)
        aprintf(32, "Invalid input polygon detected!\n");
    
    // Only take the intersection if we know that they touch.
    double baseArea;
    GEOSArea_r(ctx, base, &baseArea);

    GEOSGeometry *intersect = GEOSIntersection_r(ctx, base, overlay);
    double intersectArea;
    GEOSArea_r(ctx, intersect, &intersectArea);
    GEOSGeom_destroy_r(ctx, intersect);

    double result = intersectArea / baseArea;
    if (isnan(result))
        printf("%.4lf / %.4lf = %.4lf\n", intersectArea, baseArea, result);
    return result;
}


void getCenter(GEOSContextHandle_t ctx, GEOSGeometry *g, double *x, double *y)
{
    GEOSGeometry *center = GEOSGetCentroid_r(ctx, g);
    if (GEOSGetNumCoordinates_r(ctx, center) != 1)
    {
        aprintf(39, "Incorrect number of points in centroid\n");
        return;
    }

    GEOSCoordSeq_getXY_r(ctx, GEOSGeom_getCoordSeq_r(ctx, center), 0, x, y);
    GEOSGeom_destroy_r(ctx, center);
}


int transformer(double *x, double *y, void *userdata)
{
    double *centroid = (double *)userdata;
    // *x = (*x - centroid[0]) * 100;
//     *y = (*y - centroid[1]) * 100;
//     
     *x = (*x - centroid[0]);
    *y = (*y - centroid[1]);

    return 1;
}


GEOSGeometry *centerGeometry(GEOSContextHandle_t ctx,  GEOSGeometry *g) {
    // Clone the original geometry to ensure it remains unchanged
    GEOSGeometry* gClone = GEOSGeom_clone_r(ctx, g);

    double centroid[2];
    // Use the clone for further operations
    getCenter(ctx, gClone, &(centroid[0]), &(centroid[1]));

    GEOSGeometry *newGeo =
        GEOSGeom_transformXY_r(ctx, gClone, transformer, centroid);
    
    // Destroy the clone after use
    GEOSGeom_destroy_r(ctx, gClone);

    return newGeo; // Return the new geometry
}


// Return the MBR for a series of geometries
GEOSGeometry *minimumBoundingRectangle(GEOSContextHandle_t ctx,
                                       vector<GEOSGeometry *> geo)
{
    GEOSGeometry *boundingBox = NULL;
    GEOSGeometry *Union, *curMBR;

    for (auto cur = geo.begin(); cur != geo.end(); ++cur)
    {
        curMBR = GEOSEnvelope_r(ctx, *cur);
        if (boundingBox == NULL)
        {
            boundingBox = curMBR;
            continue;
        }

        // Take the union of the current bounding box with the current geometry
        Union = GEOSUnion_r(ctx, boundingBox, curMBR);
        GEOSGeom_destroy_r(ctx,
                           boundingBox); // Have to destroy the old box before
                                         // the new one can be created
        boundingBox = GEOSEnvelope_r(ctx, Union);

        // Free up some memory
        GEOSGeom_destroy_r(ctx, curMBR);
        GEOSGeom_destroy_r(ctx, Union);
    }

    return boundingBox;
}


// Function to calculate the union of local MBRs for a process
GEOSGeometry *calculateAggregateLocalMBR(GEOSContextHandle_t ctx, const vector<GEOSGeometry *> &localMBRs)
{   
    // Combine individual MBRs to get the local MBR for the process
    GEOSGeometry *processLocalMBR = nullptr;

    for (const auto &mbr : localMBRs)
    {       
        if (!processLocalMBR)
        {
            processLocalMBR = GEOSGeom_clone_r(ctx, mbr);
        }  
        else
        { 
            GEOSGeometry *unionResult = GEOSUnion_r(ctx, processLocalMBR, mbr);
            
            GEOSGeom_destroy_r(ctx, processLocalMBR);
            processLocalMBR = unionResult;
        }
         // Free memory for the current MBR
        GEOSGeom_destroy_r(ctx, mbr);
    }

    return processLocalMBR;
}


// Function to calculate the envelope of a geometry
GEOSGeometry *calculateEnvelope(GEOSContextHandle_t ctx, GEOSGeometry *geometry)
{
    return GEOSEnvelope_r(ctx, geometry);
}

// Return the MBR for a single geometry
GEOSGeometry *minimumBoundingRectangle(GEOSContextHandle_t ctx, GEOSGeometry *geometry)
{
    GEOSGeometry *curMBR = GEOSEnvelope_r(ctx, geometry);
    return curMBR;
}

//To create grid over MBR
vector<vector<GEOSGeometry *>> createGrid(GEOSContextHandle_t ctx,
                                          GEOSGeometry *base, int gridSize)
{
    const GEOSCoordSequence *baseCoor =
        GEOSGeom_getCoordSeq_r(ctx, GEOSGetExteriorRing_r(ctx, base));
    unsigned int length;
    GEOSCoordSeq_getSize_r(ctx, baseCoor, &length);

    // Check we got four points
    // Remember that the last one is the same as the first one so that the shape
    // is closed
    if (length != 5)
    {
        aprintf(57 + integerLength(length),
                "Returned Coordinate Sequence has unexpected dimensions! %d\n",
                length);
        return vector<vector<GEOSGeometry *>>();
    }

    double minx, maxx, x, miny, maxy, y;

    // Have to set init values separatly. While -1 would work as
    // a default for the max values, there is not sane default
    // possible for the min values.
    GEOSCoordSeq_getXY_r(ctx, baseCoor, 0, &x, &y);
    minx = x;
    maxx = x;
    miny = y;
    maxy = y;

    for (int i = 1; i < 4; i++)
    {
        // Fetch the values of the current coordinate
        GEOSCoordSeq_getXY_r(ctx, baseCoor, i, &x, &y);

        minx = x < minx ? x : minx;
        maxx = x > maxx ? x : maxx;
        miny = y < miny ? y : miny;
        maxy = y > maxy ? y : maxy;
    }

    double xRange = maxx - minx;
    double yRange = maxy - miny;

    // It feels odd that we have to do this much work to get something as simple
    // as the width and height of a shape.

    double xStep = xRange / gridSize;
    double yStep = yRange / gridSize;

    vector<vector<GEOSGeometry *>> result =
        vector<vector<GEOSGeometry *>>(gridSize);
    for (int r = 0; r < gridSize; r++)
    {
        result[r] = vector<GEOSGeometry *>(gridSize);
        for (int c = 0; c < gridSize; c++)
        {
            result[r][c] = GEOSGeom_createRectangle_r(
                ctx, minx + c * xStep, miny + r * yStep, minx + (c + 1) * xStep,
                miny + (r + 1) * yStep);
        }
    }

    return result;
}

// //to create sketch for polygon within grid
vector<double> sketch(GEOSContextHandle_t ctx,
                      vector<vector<GEOSGeometry *>> *grid, GEOSGeometry *g)

{      

    // The grid should always be square, but just in case it isn't, check the
    // second dimension too
    vector<double> result = vector<double>();

    int rowNum = grid->size();
    int colNum;
    for (int r = 0; r < rowNum; r++)
    {  //  cout<<"inside rows" <<endl;
        colNum = grid->at(r).size();
        for (int c = 0; c < colNum; c++)
        {
            //cout<<"inside columns" <<endl;
            if (grid->at(r).at(c) != NULL)
            {    //cout<<"inside grids" <<endl;
           
            //double ratio = fillRatio(ctx, grid->at(r).at(c), g);
           result.push_back(fillRatio(ctx, grid->at(r).at(c), g));
         //result.push_back(sqrt(ratio));
           // result.push_back(pow(ratio, 0.5));
                
                
            }
        }
    }

    result.shrink_to_fit();
    return result;
}





// Define the sigmoid scaling function
double sigmoidScale(double value, double k, double x0) {
    // Sigmoid scaling formula
    return 1.0 / (1.0 + exp(-k * (value - x0)));
}





vector<int> LSHHash(vector<double>* sketch, int hashLength, vector<unsigned int>* seeds) {
    // Check for invalid input
    if (!sketch || sketch->empty()) {
        throw std::invalid_argument("LSHHash received an empty sketch vector.");
    }
    if (!seeds || seeds->empty()) {
        throw std::invalid_argument("LSHHash received an empty seeds vector.");
    }

    // Step 1: Normalize the sketch vector
    vector<double> normalizedSketch(sketch->size());
    double totalSum = std::accumulate(sketch->begin(), sketch->end(), 0.0);
    if (totalSum == 0.0) {
        throw std::invalid_argument("LSHHash received a sketch vector with total sum zero.");
    }
    std::transform(sketch->begin(), sketch->end(), normalizedSketch.begin(),
                   [totalSum](double val) { return val / totalSum; });

    // Step 2: Calculate the CDF
    vector<double> cdf(normalizedSketch.size());
    std::partial_sum(normalizedSketch.begin(), normalizedSketch.end(), cdf.begin());

    // Step 3: Generate the hash
    vector<int> hash(hashLength, 0);  // Stores the retry count as the hash
    for (int i = 0; i < hashLength; i++) {
        unsigned int seed = seeds->at(i % seeds->size());  // Use seed for reproducibility
        std::mt19937 generator(seed);
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        int retryCount = 0;  // Count retries for this hash index

        while (true) {  // Retry loop for shaded region
            retryCount++;
            double randomNum = distribution(generator);

            // Map the random number to an index using the CDF
            auto it = std::lower_bound(cdf.begin(), cdf.end(), randomNum);
            int index = std::distance(cdf.begin(), it);

            // Shaded region check
            if (randomNum < std::floor(randomNum) + normalizedSketch[index]) {
                hash[i] = retryCount;  // Store retry count as the hash value
                break;  // Exit the retry loop
            }
        }
    }

    // Print the generated hash for debugging
    std::cout << "Generated Hash (Retry Counts): ";
    for (int value : hash) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return hash;
}


vector<int> IoffeHash(vector<double> *sketch, int hashLength, vector<unsigned int> *seeds)
{
    // Check for invalid input
    if (!sketch || sketch->empty()) {
        throw std::invalid_argument("IoffeHash received an empty sketch vector.");
    }
    if (!seeds || seeds->empty()) {
        throw std::invalid_argument("IoffeHash received an empty seeds vector.");
    }

    vector<int> hash(hashLength, 0);

    // Precompute random variables for each hash in hashLength
    vector<double> r(hashLength), c(hashLength), beta(hashLength);
    random_device rd;
    mt19937 gen(rd());

    // Gamma distribution (shape=2, scale=1) and uniform distribution (0, 1)
    gamma_distribution<double> gammaDist(2.0, 1.0);
    uniform_real_distribution<double> uniformDist(0.0, 1.0);

    for (int i = 0; i < hashLength; i++) {
        // Use seed for consistency
        gen.seed(seeds->at(i % seeds->size()));

        // Generate random variables
        r[i] = gammaDist(gen);
        c[i] = gammaDist(gen);
        beta[i] = uniformDist(gen);
    }

    // Iterate to generate hash values
    for (int i = 0; i < hashLength; i++) {
        double min_y = numeric_limits<double>::max();
        int min_index = 0;

        for (size_t k = 0; k < sketch->size(); k++) {
            double S_k = sketch->at(k);
            if (S_k > 0) {
                // Calculate sampling point y_k
                double y_k = exp(r[i] * (floor(log(S_k) / r[i] + beta[i]) - beta[i]));

                // Find minimum y_k
                if (y_k < min_y) {
                    min_y = y_k;
                    min_index = k;
                }
            }
        }

        // Store the index of the minimum y_k as part of the hash
        hash[i] = min_index;
    }

    return hash;
}


// Function to check if a point is inside a polygon using ray-casting
bool isPointInsidePolygon(const Point& p, const vector<Point>& polygon) {
    int intersections = 0;
    size_t n = polygon.size();
    for (size_t i = 0; i < n; ++i) {
        Point v1 = polygon[i];
        Point v2 = polygon[(i + 1) % n]; // Wrap around to the first vertex

        // Check if the ray crosses the edge
        if ((v1.y > p.y) != (v2.y > p.y)) {
            double xIntersection = v1.x + (p.y - v1.y) * (v2.x - v1.x) / (v2.y - v1.y);
            if (p.x < xIntersection) {
                intersections++;
            }
        }
    }
    return intersections % 2 == 1; // Odd intersections mean inside
}



// //based on global MBR not radial
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);  // Result hash vector
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute the local MBR of the polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
// 
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
// 
//     // Step 3: Set up uniform distributions over the global MBR
//     uniform_real_distribution<double> xDist(bbox.getXMin(), bbox.getXMax());
//     uniform_real_distribution<double> yDist(bbox.getYMin(), bbox.getYMax());
// 
//     // Step 4: Precompute inverse grid cell size
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; ++i) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             ++attempts;
// 
//             // Step 5: Generate random point uniformly in global MBR
//             double dartX = xDist(generator);
//             double dartY = yDist(generator);
// 
//             // Step 6: Quickly discard if outside local polygon's MBR
//             if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax)
//                 continue;
// 
//             // Step 7: Map point to grid cell
//             int gridX = static_cast<int>((dartX - bbox.getXMin()) * invCellSizeX);
//             int gridY = static_cast<int>((dartY - bbox.getYMin()) * invCellSizeY);
// 
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY)
//                 continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
// 
//             // Step 8: Reject if cell is not occupied
//             if (occupiedCells.find(cellIndex) == occupiedCells.end())
//                 continue;
// 
//             // Step 9: Final point-in-polygon test
//             if (isPointInsidePolygon({dartX, dartY}, polygon)) {
//                 hash[i] = attempts;
//                 break;
//             }
//         }
//     }
// 
//     return hash;
// } 

// //reduced mbr space by 10%
// vector<int> hash2d(
//     const vector<Point>& polygon,
//     int hashLength,
//     const BoundingBox& bbox,
//     const vector<vector<GEOSGeometry*>>& grid,
//     const std::unordered_set<int>& occupiedCells,
//     const vector<unsigned int>& seeds,
//     GEOSContextHandle_t ctx
// ) {
//     int gridSizeX = (int)grid.size();
//     int gridSizeY = (int)grid[0].size();
// 
//     // 1) Compute the polygon's local MBR
//     double polyXMin =  std::numeric_limits<double>::infinity();
//     double polyXMax = -std::numeric_limits<double>::infinity();
//     double polyYMin =  std::numeric_limits<double>::infinity();
//     double polyYMax = -std::numeric_limits<double>::infinity();
//     for (auto& p : polygon) {
//         polyXMin = std::min(polyXMin, p.x);
//         polyXMax = std::max(polyXMax, p.x);
//         polyYMin = std::min(polyYMin, p.y);
//         polyYMax = std::max(polyYMax, p.y);
//     }
// 
//     // 2) Build and prepare the GEOS polygon
//     GEOSCoordSequence* cs = GEOSCoordSeq_create_r(ctx, polygon.size()+1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, cs, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, cs, i, polygon[i].y);
//     }
//     // close the ring
//     GEOSCoordSeq_setX_r(ctx, cs, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, cs, polygon.size(), polygon[0].y);
// 
//     GEOSGeometry* ring    = GEOSGeom_createLinearRing_r(ctx, cs);
//     GEOSGeometry* gPoly   = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepPoly = GEOSPrepare_r(ctx, gPoly);
// 
//     // 3) Compute a 10%-shrunken global MBR for sampling
//     double X0 = bbox.getXMin();
//     double X1 = bbox.getXMax();
//     double Y0 = bbox.getYMin();
//     double Y1 = bbox.getYMax();
//     double dx = X1 - X0;
//     double dy = Y1 - Y0;
//     // shrink 10% on each side
//     double sampleX0 = X0 + 0.1 * dx;
//     double sampleX1 = X1 - 0.1 * dx;
//     double sampleY0 = Y0 + 0.1 * dy;
//     double sampleY1 = Y1 - 0.1 * dy;
// 
//     // 4) Prepare RNGs and distributions over the shrunken MBR
//     vector<std::mt19937> gens;
//     for (auto s : seeds) gens.emplace_back(s);
//     std::uniform_real_distribution<double> xDist(sampleX0, sampleX1);
//     std::uniform_real_distribution<double> yDist(sampleY0, sampleY1);
// 
//     // 5) Precompute inverse cell size for mapping (x,y) → cell index
//     double invCellX = gridSizeX / (X1 - X0);
//     double invCellY = gridSizeY / (Y1 - Y0);
// 
//     // 6) Main hashing loop
//     vector<int> hash(hashLength, 0);
//     for (int i = 0; i < hashLength; ++i) {
//         auto& gen = gens[i % gens.size()];
//         int attempts = 0;
// 
//         while (true) {
//             ++attempts;
//             double x = xDist(gen), y = yDist(gen);
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, x, y);
// 
//             // a) Quickly reject if outside the polygon's local MBR
//             if (x < polyXMin || x > polyXMax ||
//                 y < polyYMin || y > polyYMax) {
//                 GEOSGeom_destroy_r(ctx, pt);
//                 continue;
//             }
// 
//             // b) Test real polygon containment
//             if (GEOSPreparedContains_r(ctx, prepPoly, pt)) {
//                 // Map to grid cell in the original global MBR
//                 int gx = std::min(std::max(int((x - X0) * invCellX), 0), gridSizeX-1);
//                 int gy = std::min(std::max(int((y - Y0) * invCellY), 0), gridSizeY-1);
//                 int cellIdx = gy * gridSizeX + gx;
// 
//                 // Only accept if that cell is occupied
//                 if (occupiedCells.count(cellIdx)) {
//                     hash[i] = attempts;
//                     GEOSGeom_destroy_r(ctx, pt);
//                     break;
//                 }
//             }
// 
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     // 7) Cleanup GEOS objects
//     GEOSPreparedGeom_destroy_r(ctx, prepPoly);
//     GEOSGeom_destroy_r(ctx, gPoly);
// 
//     return hash;
// }


// // based on global MBR not radial
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                     GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);  // Result hash vector
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute the local MBR of the polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
// 
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     
//     
//     // Step 4: Convert polygon to GEOSGeometry and prepare it
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     // Step 3: Set up uniform distributions over the global MBR
//     uniform_real_distribution<double> xDist(bbox.getXMin(), bbox.getXMax());
//     uniform_real_distribution<double> yDist(bbox.getYMin(), bbox.getYMax());
// 
//     // Step 4: Precompute inverse grid cell size
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; ++i) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             ++attempts;
// 
//             // Step 5: Generate random point uniformly in global MBR
//             double dartX = xDist(generator);
//             double dartY = yDist(generator);
// 
//             // Step 6: Quickly discard if outside local polygon's MBR
//             if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax)
//                 continue;
// 
//             // Step 7: Map point to grid cell
//             int gridX = static_cast<int>((dartX - bbox.getXMin()) * invCellSizeX);
//             int gridY = static_cast<int>((dartY - bbox.getYMin()) * invCellSizeY);
// 
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY)
//                 continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
// 
//             // Step 8: Reject if cell is not occupied
//             if (occupiedCells.find(cellIndex) == occupiedCells.end())
//                 continue;
// 
//              GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     return hash;
// } 

// based on global MBR not radial
vector<int> hash2d(const vector<Point>& polygon, int hashLength,
                   const BoundingBox& bbox,
                   const vector<vector<GEOSGeometry*>>& grid,
                   const unordered_set<int>& occupiedCells,
                   const vector<unsigned int>& seeds,
                    GEOSContextHandle_t ctx) {

    vector<int> hash(hashLength, 0);  // Result hash vector


    // Step 1: Compute the local MBR of the polygon
    double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
    double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();

    for (const auto& point : polygon) {
        polyXMin = min(polyXMin, point.x);
        polyXMax = max(polyXMax, point.x);
        polyYMin = min(polyYMin, point.y);
        polyYMax = max(polyYMax, point.y);
    }

    // Step 2: Initialize RNGs
    vector<std::mt19937> generators;
    for (unsigned int seed : seeds) {
        generators.emplace_back(seed);
    }
    
    
    // Step 4: Convert polygon to GEOSGeometry and prepare it
    GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
    for (size_t i = 0; i < polygon.size(); ++i) {
        GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
        GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
    }
    GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
    GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
    GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
    GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
    const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);

    // Step 3: Set up uniform distributions over the global MBR
    uniform_real_distribution<double> xDist(bbox.getXMin(), bbox.getXMax());
    uniform_real_distribution<double> yDist(bbox.getYMin(), bbox.getYMax());


    for (int i = 0; i < hashLength; ++i) {
        int attempts = 0;
        std::mt19937& generator = generators[i % seeds.size()];

        while (true) {
            ++attempts;

            // Step 5: Generate random point uniformly in global MBR
            double dartX = xDist(generator);
            double dartY = yDist(generator);

            // Step 6: Quickly discard if outside local polygon's MBR
            if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax)
                continue;


             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
            if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
                hash[i] = attempts;
                GEOSGeom_destroy_r(ctx, pt);
                break;
            }
            GEOSGeom_destroy_r(ctx, pt);
        }
    }

    return hash;
} 


// //diamond
// vector<int> hash2d(
//     const vector<Point>& polygon,
//     int hashLength,
//     const BoundingBox& bbox,
//     const vector<vector<GEOSGeometry*>>& grid,
//     const unordered_set<int>& occupiedCells,
//     const vector<unsigned int>& seeds,
//     GEOSContextHandle_t ctx
// ) {
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // 1) Compute local polygon MBR
//     double polyXMin =  std::numeric_limits<double>::max();
//     double polyXMax = -std::numeric_limits<double>::max();
//     double polyYMin =  std::numeric_limits<double>::max();
//     double polyYMax = -std::numeric_limits<double>::max();
//     for (auto& p : polygon) {
//         polyXMin = std::min(polyXMin, p.x);
//         polyXMax = std::max(polyXMax, p.x);
//         polyYMin = std::min(polyYMin, p.y);
//         polyYMax = std::max(polyYMax, p.y);
//     }
// 
//     // 2) Build GEOS polygon + prepare
//     GEOSCoordSequence* seq = GEOSCoordSeq_create_r(ctx, polygon.size()+1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, seq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, seq, i, polygon[i].y);
//     }
//     // close ring
//     GEOSCoordSeq_setX_r(ctx, seq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, seq, polygon.size(), polygon[0].y);
// 
//     GEOSGeometry* ring       = GEOSGeom_createLinearRing_r(ctx, seq);
//     GEOSGeometry* geosPoly    = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepPoly = GEOSPrepare_r(ctx, geosPoly);
// 
//     // 3) Build midpoint‐diamond inside global MBR
//     double X0 = bbox.getXMin(), X1 = bbox.getXMax();
//     double Y0 = bbox.getYMin(), Y1 = bbox.getYMax();
//     double mX = 0.5*(X0+X1), mY = 0.5*(Y0+Y1);
// 
//     // diamond coords (midpoints of each edge)
//     vector<Point> diamond = {
//         { X0, mY },
//         { mX, Y1 },
//         { X1, mY },
//         { mX, Y0 },
//         { X0, mY }  // close ring
//     };
//     GEOSCoordSequence* dseq = GEOSCoordSeq_create_r(ctx, diamond.size(), 2);
//     for (size_t i = 0; i < diamond.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, dseq, i, diamond[i].x);
//         GEOSCoordSeq_setY_r(ctx, dseq, i, diamond[i].y);
//     }
//     GEOSGeometry* dring        = GEOSGeom_createLinearRing_r(ctx, dseq);
//     GEOSGeometry* midRectPoly  = GEOSGeom_createPolygon_r(ctx, dring, nullptr, 0);
//     const GEOSPreparedGeometry* prepMid = GEOSPrepare_r(ctx, midRectPoly);
// 
//     // 4) RNGs & distributions over the **global** MBR
//     vector<mt19937> gens;
//     for (auto s : seeds) gens.emplace_back(s);
//     uniform_real_distribution<double> xDist(X0, X1);
//     uniform_real_distribution<double> yDist(Y0, Y1);
// 
//     // 5) Inverse cell sizes for mapping to grid
//     double invCellX = gridSizeX / (X1 - X0);
//     double invCellY = gridSizeY / (Y1 - Y0);
// 
//     vector<int> hash(hashLength);
//     for (int i = 0; i < hashLength; ++i) {
//         auto& gen = gens[i % gens.size()];
//         int attempts = 0;
// 
//         while (true) {
//             ++attempts;
//             double x = xDist(gen), y = yDist(gen);
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, x, y);
// 
//             // First reject if outside the midpoint‐diamond
//             if (!GEOSPreparedContains_r(ctx, prepMid, pt)) {
//                 GEOSGeom_destroy_r(ctx, pt);
//                 continue;
//             }
// 
//             // Then test actual polygon
//             if (GEOSPreparedContains_r(ctx, prepPoly, pt)) {
//                 // map to grid cell
//                 int gx = std::min(std::max(int((x - X0)*invCellX), 0), gridSizeX-1);
//                 int gy = std::min(std::max(int((y - Y0)*invCellY), 0), gridSizeY-1);
//                 int cellIdx = gy * gridSizeX + gx;
//                 hash[i] = cellIdx;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
// 
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     // 6) Cleanup
//     GEOSPreparedGeom_destroy_r(ctx, prepPoly);
//     GEOSGeom_destroy_r(ctx, geosPoly);
//     GEOSPreparedGeom_destroy_r(ctx, prepMid);
//     GEOSGeom_destroy_r(ctx, midRectPoly);
// 
//     return hash;
// }
// 
// 
// //traingle
// 
// vector<int> hash2d(
//     const vector<Point>& polygon,
//     int hashLength,
//     const BoundingBox& bbox,
//     int gridSizeX, int gridSizeY,
//     const vector<unsigned int>& seeds,
//     GEOSContextHandle_t ctx
// ) {
//     // 1) Build GEOS polygon + prepare
//     GEOSCoordSequence* cs = GEOSCoordSeq_create_r(ctx, polygon.size()+1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, cs, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, cs, i, polygon[i].y);
//     }
//     // close the ring
//     GEOSCoordSeq_setX_r(ctx, cs, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, cs, polygon.size(), polygon[0].y);
// 
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, cs);
//     GEOSGeometry* geoPoly = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     // (optional) prepared for fast contains tests:
//     const GEOSPreparedGeometry* prepPoly = GEOSPrepare_r(ctx, geoPoly);
// 
//     // 2) Triangulate the polygon
//     vector<GEOSGeometry*> triangles = GEOSPolygonTriangulate_r(ctx, geoPoly);
//     size_t nTri = triangles.size();
// 
//     // 3) Compute areas & cumulative distribution
//     vector<double> areas(nTri);
//     double totalArea = 0;
//     for (size_t t = 0; t < nTri; ++t) {
//         areas[t] = GEOSArea_r(ctx, triangles[t]);
//         totalArea += areas[t];
//     }
//     vector<double> cumArea(nTri);
//     double c = 0;
//     for (size_t t = 0; t < nTri; ++t) {
//         c += areas[t];
//         cumArea[t] = c / totalArea;  // will end at 1.0
//     }
// 
//     // 4) Prepare RNGs & distributions
//     vector<mt19937> gens;
//     for (auto s : seeds) gens.emplace_back(s);
//     uniform_real_distribution<double> u01(0.0, 1.0);
// 
//     double width  = bbox.getXMax() - bbox.getXMin();
//     double height = bbox.getYMax() - bbox.getYMin();
//     double invCellX = gridSizeX / width;
//     double invCellY = gridSizeY / height;
// 
//     // 5) For each hash dimension, sample one point & record its cell
//     vector<int> hash(hashLength);
//     for (int i = 0; i < hashLength; ++i) {
//         mt19937& gen = gens[i % gens.size()];
// 
//         // a) pick triangle by area weight
//         double u = u01(gen);
//         size_t triIdx = std::distance(
//             cumArea.begin(),
//             upper_bound(cumArea.begin(), cumArea.end(), u)
//         );
// 
//         // b) sample inside that triangle via barycentric coords
//         double r1 = u01(gen), r2 = u01(gen);
//         double srt = std::sqrt(r1);
//         double b0 = 1 - srt;
//         double b1 = srt * (1 - r2);
//         double b2 = srt * r2;
// 
//         // get the 3 vertices of triangles[triIdx]
//         GEOSCoordSequence* tcs = GEOSGeom_getCoordSeq_r(ctx, triangles[triIdx]);
//         double x0,y0,x1,y1,x2,y2;
//         GEOSCoordSeq_getX_r(ctx, tcs, 0, &x0);  GEOSCoordSeq_getY_r(ctx, tcs, 0, &y0);
//         GEOSCoordSeq_getX_r(ctx, tcs, 1, &x1);  GEOSCoordSeq_getY_r(ctx, tcs, 1, &y1);
//         GEOSCoordSeq_getX_r(ctx, tcs, 2, &x2);  GEOSCoordSeq_getY_r(ctx, tcs, 2, &y2);
// 
//         double px = b0*x0 + b1*x1 + b2*x2;
//         double py = b0*y0 + b1*y1 + b2*y2;
// 
//         // c) map (px,py) to grid cell index
//         int gx = int((px - bbox.getXMin()) * invCellX);
//         int gy = int((py - bbox.getYMin()) * invCellY);
//         // clamp just in case of fp error
//         gx = std::min(std::max(gx, 0), gridSizeX - 1);
//         gy = std::min(std::max(gy, 0), gridSizeY - 1);
// 
//         hash[i] = gy * gridSizeX + gx;
//     }
// 
//     // 6) Cleanup GEOS objects
//     GEOSPreparedGeom_destroy_r(ctx, prepPoly);
//     GEOSGeom_destroy_r(ctx, geoPoly);
//     for (auto t : triangles) GEOSGeom_destroy_r(ctx, t);
// 
//     return hash;
// }

// vector<int> hash2d(const vector<Point>& polygon, int hashLength, vector<unsigned int>* seeds,
//                    double xMin, double xMax, double yMin, double yMax) {
//     vector<int> hash(hashLength, 0); // Result hash vector
// 
//     // Compute the local MBR of the polygon (tight bounding box)
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
// 
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Use the global MBR for sampling
//     uniform_real_distribution<double> xDist(xMin, xMax);
//     uniform_real_distribution<double> yDist(yMin, yMax);
// 
//     // Generate hashes
//     for (int i = 0; i < hashLength; i++) {
//         minstd_rand generator(seeds->at(i % seeds->size()));
//         int attempts = 0;
// 
//         while (true) {
//             attempts++; // Count number of tries
//             Point dart = {xDist(generator), yDist(generator)}; // Random dart
// 
//             // Check if the point is inside the polygon's local MBR first
//             if (dart.x < polyXMin || dart.x > polyXMax || dart.y < polyYMin || dart.y > polyYMax) {
//                 continue; // Skip and try again
//             }
// 
//             // Check if the point is inside the polygon only if it's inside the local MBR
//             if (isPointInsidePolygon(dart, polygon)) {
//                 hash[i] = attempts; // Store attempts needed for the dart to hit
//                 break; // Stop once the dart lands inside
//             }
//         }
//     }
// 
//     return hash;
// }




BoundingBox extractBoundingBox(GEOSContextHandle_t ctx, GEOSGeometry* envelope) {
    double xMin, yMin, xMax, yMax;

    // Check if the envelope is valid and retrieve bounds
    if (!GEOSGeom_getExtent_r(ctx, envelope, &xMin, &yMin, &xMax, &yMax)) {
        throw std::runtime_error("Failed to retrieve bounds for the global MBR.");
    }

    // Return a BoundingBox object
    return BoundingBox(xMin, yMin, xMax, yMax);
}



// Helper structure for grid points
struct GridPoint {
    int x, y;
};




unordered_set<int> precomputeOccupiedCells(const vector<double>& sketch, int gridSizeX, int gridSizeY) {
    unordered_set<int> occupiedCells;
    for (int i = 0; i < sketch.size(); i++) {
        if (sketch[i] > 0) occupiedCells.insert(i);
    }
    return occupiedCells;
}



// // uniform radial hash
// vector<int> hash2d(const vector<Point>& polygon, int hashLength, 
//                    const BoundingBox& bbox, 
//                    const vector<vector<GEOSGeometry*>>& grid, 
//                    const unordered_set<int>& occupiedCells,  // Precomputed occupied grid cells
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// // ) {  // Pass seeds as a parameter
// 
//     vector<int> hash(hashLength, 0);  // Result hash vector
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
// 
//      // Step 1: Compute the local MBR of the polygon
//      // double polyMBRStart = MPI_Wtime();  // Start timing
// 
//     // Step 1: Compute the local MBR of the polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
// 
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// //     
// //     double polyMBREnd = MPI_Wtime();  // End timing
// // double polyMBRTime = polyMBREnd - polyMBRStart;
// 
// // Only print from rank 0, or if you prefer, remove rank check if not available in this scope
// //std::cout << "[hash2d] Polygon MBR time: " << polyMBRTime << " seconds\n";
// 
// 
//     // Step 2: Compute the **center** of the global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double maxRadius = max(bbox.getXMax() - centerX, bbox.getYMax() - centerY);  // Furthest point from center
// 
//     // Step 3: Initialize random generators
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed); // Create RNGs with fixed seeds
//     }
// 
//      //  // Step 4: Distributions for **radius and angle** (r, θ)
// //        uniform_real_distribution<double> thetaDist(0, 2 * M_PI);  // Random angle (0 to 360 degrees)
// //        uniform_real_distribution<double> radiusDist(0, maxRadius);  // Expanding radius
// 
//        // Step 4: Distributions for **radius and angle** (r, θ)
//        uniform_real_distribution<double> thetaDist(0, 2 * M_PI);  // Random angle (0 to 360 degrees)
//        uniform_real_distribution<double> uniformDist(0, 1);  // Corrected: Added uniform distribution for bias correction
// 
//     // Step 5: Compute grid cell size to avoid repeated calculations
//     // double cellSizeX = (bbox.getXMax() - bbox.getXMin()) / gridSizeX;
//     //double cellSizeY = (bbox.getYMax() - bbox.getYMin()) / gridSizeY;
//       double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//       double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];  // Use indexed seed
// 
//         while (true) {
//             attempts++;  // Count attempts
// 
//              // **Step 6: Generate a random point with bias correction**
//             double u = uniformDist(generator);   // Uniform sample in [0,1]
//             double r = maxRadius * sqrt(u);   // Corrected radius distribution
//             double theta = thetaDist(generator);  // Get random angle
//  
// 
//             // Convert (r, θ) → (x, y)
//             double dartX = centerX + r * cos(theta);
//             double dartY = centerY + r * sin(theta);
// 
//             // **Step 7: First check if the point is inside the polygon's local MBR**
//             // if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax) {
// //                 continue;  // Skip and retry
// //             }
//              bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                    dartY >= polyYMin && dartY <= polyYMax);
//              if (!withinPolyMBR) continue;
// 
//             // **Step 8: Map the random point to a grid cell (Only if inside local MBR)**
// //             int gridX = (dartX - bbox.getXMin()) / cellSizeX;
// //             int gridY = (dartY - bbox.getYMin()) / cellSizeY;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             
//             // Ensure the grid indices are within bounds
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) {
//                 continue;  // Retry if out of bounds
//             }
// 
//             // **Step 9: Compute the linear index of the grid cell**
//             int cellIndex = gridY * gridSizeX + gridX;
// 
//             // **Step 10: Skip empty cells (Use precomputed occupied cells)**
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             // **Step 11: Perform PnP test only if necessary**
//             if (isPointInsidePolygon({dartX, dartY}, polygon)) {
//                 hash[i] = attempts;  // Store attempts needed for the dart to land inside
//                 break;  // Stop once the dart lands inside
//             }
//         }
//     }
// 
//     return hash;
// }

// 
// //optimized hash2d
// // Optimized hash2d using GEOSPreparedContains
// vector<int> hash2d(const vector<Point>& polygon, int hashLength, 
//                    const BoundingBox& bbox, 
//                    const vector<vector<GEOSGeometry*>>& grid, 
//                    const unordered_set<int>& occupiedCells, 
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute local MBR of polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Compute center of global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double maxRadius = max(bbox.getXMax() - centerX, bbox.getYMax() - centerY);
// 
//     // Step 3: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
//     uniform_real_distribution<double> uniformDist(0, 1);
// 
//     // Step 4: Convert polygon to GEOSGeometry and prepare it
//     GEOSGeometry* geosPolygon = createGEOSPolygonFromPoints(ctx, polygon);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             attempts++;
//             double u = uniformDist(generator);
//             double r = maxRadius * sqrt(u);
//             double theta = thetaDist(generator);
// 
//             double dartX = centerX + r * cos(theta);
//             double dartY = centerY + r * sin(theta);
// 
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }

// // Optimized hash2d using GEOSPreparedContains
// vector<int> hash2d(const vector<Point>& polygon, int hashLength, 
//                    const BoundingBox& bbox, 
//                    const vector<vector<GEOSGeometry*>>& grid, 
//                    const unordered_set<int>& occupiedCells, 
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
//                    
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute local MBR of polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Compute center of global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double maxRadius = max(bbox.getXMax() - centerX, bbox.getYMax() - centerY);
// // 
// // //✅ Diagonal radius from center to corner
// //     double centerX = (bbox.getXMax() - bbox.getXMin()) / 2.0;
// //     double centerY = (bbox.getYMax() - bbox.getYMin()) / 2.0;
// //     double maxRadius = sqrt(centerX * centerX + centerY * centerY);
// // 
// 
//     // Step 3: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
//     uniform_real_distribution<double> uniformDist(0, 1);
// 
//     // Step 4: Convert polygon to GEOSGeometry and prepare it
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             attempts++;
//             double u = uniformDist(generator);
//             double r = maxRadius * sqrt(u);
//             double theta = thetaDist(generator);
// 
//             double dartX = centerX + r * cos(theta);
//             double dartY = centerY + r * sin(theta);
// 
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }


// // ✅ Disk → replaced with ellipse inside MBR
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute local MBR of polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // ✅ Step 2: ellipse inside global MBR
//     double width = bbox.getXMax() - bbox.getXMin();
//     double height = bbox.getYMax() - bbox.getYMin();
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double a = width / 2.0;   // semi-axis X
//     double b = height / 2.0;  // semi-axis Y
// 
//     // Step 3: RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
//     uniform_real_distribution<double> uniformDist(0, 1);
// 
//     // Step 4: Convert polygon to GEOS
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / width;
//     double invCellSizeY = gridSizeY / height;
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             attempts++;
//             double u = uniformDist(generator);
//             double sqrtU = sqrt(u);
//             double theta = thetaDist(generator);
// 
//             // ✅ sample ellipse
//             double dartX = centerX + a * sqrtU * cos(theta);
//             double dartY = centerY + b * sqrtU * sin(theta);
// 
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }

// // ✅ Sampling inside MBR in "rectangular polar" coordinates
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Polygon MBR
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // ✅ Step 2: MBR center
//     double width = bbox.getXMax() - bbox.getXMin();
//     double height = bbox.getYMax() - bbox.getYMin();
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double a = width / 2.0;
//     double b = height / 2.0;
// 
//     // Step 3: RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> uniformX(-1, 1);  // scaled
//     uniform_real_distribution<double> uniformY(-1, 1);
// 
//     // Step 4: Convert polygon to GEOS
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / width;
//     double invCellSizeY = gridSizeY / height;
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             attempts++;
// 
//             // ✅ Uniformly sample in [-1, 1] and scale to MBR
//             double ux = uniformX(generator);
//             double uy = uniformY(generator);
// 
//             double dartX = centerX + ux * a;
//             double dartY = centerY + uy * b;
// 
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }

// // 🏆 GLOBAL lookup table: index → list of (u, theta, cos, sin) pairs
// static unordered_map<int, vector<tuple<double, double, double, double>>> lookupTable;
// 
// // Optimized hash2d with per-index cached (u, theta, cos, sin)
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute polygon MBR
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Global MBR center
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double maxRadius = max(bbox.getXMax() - centerX, bbox.getYMax() - centerY);
// 
// // //✅ Diagonal radius from center to corner
// //     double centerX = (bbox.getXMax() - bbox.getXMin()) / 2.0;
// //     double centerY = (bbox.getYMax() - bbox.getYMin()) / 2.0;
// //     double maxRadius = sqrt(centerX * centerX + centerY * centerY);
// 
// 
//     // Step 3: RNGs
//     vector<mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
//     uniform_real_distribution<double> uniformDist(0, 1);
// 
//     // Step 4: Convert polygon to GEOS
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             // Check if we already cached this (index, attempt)
//             if (lookupTable[i].size() <= attempts) {
//                 // If not cached → generate new (u, theta) and store
//                 double u = uniformDist(generator);
//                 double theta = thetaDist(generator);
//                 lookupTable[i].emplace_back(u, theta, cos(theta), sin(theta));
//             }
// 
//             // Retrieve cached values
//             const auto& [u, theta, cosTheta, sinTheta] = lookupTable[i][attempts];
//             double r = maxRadius * sqrt(u);
//             double dartX = centerX + r * cosTheta;
//             double dartY = centerY + r * sinTheta;
// 
//             attempts++;
// 
//             // Optional bounding box check
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }
//using W and L and precomputed trigonometric tables// 
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute local MBR of polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Compute center of global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
// 
//     // ✅ Diagonal radius from center to corner
//     double halfWidth = (bbox.getXMax() - bbox.getXMin()) / 2.0;
//     double halfHeight = (bbox.getYMax() - bbox.getYMin()) / 2.0;
//     double maxRadius = sqrt(halfWidth * halfWidth + halfHeight * halfHeight);
// 
//     // Step 3: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> uniformDist(0, 1);
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
// 
//     // ✅ Precompute cos/sin per hash dimension (instead of table)
//     vector<pair<double, double>> precomputedAngles(hashLength);
//     for (int i = 0; i < hashLength; ++i) {
//         std::mt19937& generator = generators[i % seeds.size()];
//         double theta = thetaDist(generator);
//         precomputedAngles[i] = { cos(theta), sin(theta) };
//     }
// 
//     // Step 4: Convert polygon to GEOS geometry
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
//         const auto& angle = precomputedAngles[i];  // cached cos/sin
// 
//         while (true) {
//             attempts++;
//             double u = uniformDist(generator);
//             double r = maxRadius * sqrt(u);  // radial sampling
// 
//             double dartX = centerX + r * angle.first;
//             double dartY = centerY + r * angle.second;
// 
//             // Optional: quick bounding check
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }
// // Optimized hash2d using GEOSPreparedContains using both L and W not max
// vector<int> hash2d(const vector<Point>& polygon, int hashLength, 
//                    const BoundingBox& bbox, 
//                    const vector<vector<GEOSGeometry*>>& grid, 
//                    const unordered_set<int>& occupiedCells, 
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
//                    
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute local MBR of polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMin, point.y);
//     }
// 
//     // Step 2: Compute center of global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
// 
//     // ✅ Updated: compute radius using diagonal (distance from center to corner)
//     double halfWidth = (bbox.getXMax() - bbox.getXMin()) / 2.0;
//     double halfHeight = (bbox.getYMax() - bbox.getYMin()) / 2.0;
//     double maxRadius = sqrt(halfWidth * halfWidth + halfHeight * halfHeight);
// 
//     // Step 3: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
//     uniform_real_distribution<double> uniformDist(0, 1);
// 
//     // Step 4: Convert polygon to GEOSGeometry and prepare it
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             attempts++;
//             double u = uniformDist(generator);
//             double r = maxRadius * sqrt(u);  // radial sampling
//             double theta = thetaDist(generator);
// 
//             double dartX = centerX + r * cos(theta);
//             double dartY = centerY + r * sin(theta);
// 
//             // Optional: quick bounding check with polygon's local MBR
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }
// 

// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute local MBR of polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Compute center of global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
// 
//     // ✅ Diagonal radius from center to corner
//     double halfWidth = (bbox.getXMax() - bbox.getXMin()) / 2.0;
//     double halfHeight = (bbox.getYMax() - bbox.getYMin()) / 2.0;
//     double maxRadius = sqrt(halfWidth * halfWidth + halfHeight * halfHeight);
// 
//     // Step 3: Static precomputed cos/sin table
//     static vector<double> precomputedCos;
//     static vector<double> precomputedSin;
//     const int NUM_ANGLES = 1024;
//     if (precomputedCos.empty()) {
//         precomputedCos.resize(NUM_ANGLES);
//         precomputedSin.resize(NUM_ANGLES);
//         for (int j = 0; j < NUM_ANGLES; ++j) {
//             double angle = 2.0 * M_PI * j / NUM_ANGLES;
//             precomputedCos[j] = cos(angle);
//             precomputedSin[j] = sin(angle);
//         }
//     }
// 
//     // Step 4: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
//     uniform_real_distribution<double> uniformDist(0, 1);
//     uniform_int_distribution<int> angleIndexDist(0, NUM_ANGLES - 1);
// 
//     // Step 5: Convert polygon to GEOS geometry
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; i++) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             attempts++;
//             double u = uniformDist(generator);
//             double r = maxRadius * sqrt(u);  // radial sampling
//             int thetaIdx = angleIndexDist(generator);  // random precomputed angle
// 
//             double dartX = centerX + r * precomputedCos[thetaIdx];
//             double dartY = centerY + r * precomputedSin[thetaIdx];
// 
//             // Optional: quick bounding check
//             bool withinPolyMBR = (dartX >= polyXMin && dartX <= polyXMax &&
//                                   dartY >= polyYMin && dartY <= polyYMax);
//             if (!withinPolyMBR) continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY) continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     GEOSPreparedGeom_destroy_r(ctx, prepGeom);
//     GEOSGeom_destroy_r(ctx, geosPolygon);
//     return hash;
// }



// //starts from the center
// vector<int> hash2d(const vector<Point>& polygon, int hashLength, 
//                    const BoundingBox& bbox, 
//                    const vector<vector<GEOSGeometry*>>& grid, 
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds) {
// 
//     vector<int> hash(hashLength, 0);
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute the polygon's local MBR
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Global center and max radius of bounding disk
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double maxRadius = max(bbox.getXMax() - centerX, bbox.getYMax() - centerY);
// 
//     // Step 3: Random number generators (one per seed)
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
// 
//     uniform_real_distribution<double> thetaDist(0, 2 * M_PI);
// 
//     // Step 4: Precompute inverse grid cell sizes
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     // Step 5: Dart throwing from center outward (deterministic radii)
//     for (int i = 0; i < hashLength; ++i) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         const int maxSteps = 10000; // Max concentric layers
//         double radiusStep = maxRadius / maxSteps;
// 
//         for (int step = 1; step <= maxSteps; ++step) {
//             ++attempts;
// 
//             double r = step * radiusStep;
//             double theta = thetaDist(generator);
// 
//             double dartX = centerX + r * cos(theta);
//             double dartY = centerY + r * sin(theta);
// 
//             // Skip if outside local MBR
//             if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax)
//                 continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
// 
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY)
//                 continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
//             if (occupiedCells.find(cellIndex) == occupiedCells.end()) continue;
// 
//             if (isPointInsidePolygon({dartX, dartY}, polygon)) {
//                 hash[i] = attempts;
//                 break;
//             }
//         }
// 
//         // If maxSteps exhausted, hash[i] stays 0 (or you can assign a sentinel value like maxSteps + 1)
//     }
// 
//     return hash;
// }



// // based on global MBR not radial optimized
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds,
//                    GEOSContextHandle_t ctx) {
// 
//     vector<int> hash(hashLength, 0);  // Result hash vector
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Compute the local MBR of the polygon
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
// 
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
// 
//     // Step 3: Set up uniform distributions over the global MBR
//     uniform_real_distribution<double> xDist(bbox.getXMin(), bbox.getXMax());
//     uniform_real_distribution<double> yDist(bbox.getYMin(), bbox.getYMax());
//     
//    // Step 4: Convert polygon to GEOS
//     GEOSCoordSequence* coordSeq = GEOSCoordSeq_create_r(ctx, polygon.size() + 1, 2);
//     for (size_t i = 0; i < polygon.size(); ++i) {
//         GEOSCoordSeq_setX_r(ctx, coordSeq, i, polygon[i].x);
//         GEOSCoordSeq_setY_r(ctx, coordSeq, i, polygon[i].y);
//     }
//     
//     GEOSCoordSeq_setX_r(ctx, coordSeq, polygon.size(), polygon[0].x);
//     GEOSCoordSeq_setY_r(ctx, coordSeq, polygon.size(), polygon[0].y);
//     GEOSGeometry* ring = GEOSGeom_createLinearRing_r(ctx, coordSeq);
//     GEOSGeometry* geosPolygon = GEOSGeom_createPolygon_r(ctx, ring, nullptr, 0);
//     const GEOSPreparedGeometry* prepGeom = GEOSPrepare_r(ctx, geosPolygon);
// 
//     // Step 4: Precompute inverse grid cell size
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     for (int i = 0; i < hashLength; ++i) {
//         int attempts = 0;
//         std::mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             ++attempts;
// 
//             // Step 5: Generate random point uniformly in global MBR
//             double dartX = xDist(generator);
//             double dartY = yDist(generator);
// 
//             // Step 6: Quickly discard if outside local polygon's MBR
//             if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax)
//                 continue;
// 
//             // Step 7: Map point to grid cell
//             int gridX = static_cast<int>((dartX - bbox.getXMin()) * invCellSizeX);
//             int gridY = static_cast<int>((dartY - bbox.getYMin()) * invCellSizeY);
// 
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY)
//                 continue;
// 
//             // int cellIndex = gridY * gridSizeX + gridX;
// // 
// //             // Step 8: Reject if cell is not occupied
// //             if (occupiedCells.find(cellIndex) == occupiedCells.end())
// //                 continue;
// 
//             // Step 9: Final point-in-polygon test
//             GEOSGeometry* pt = GEOSGeom_createPointFromXY_r(ctx, dartX, dartY);
//             if (GEOSPreparedContains_r(ctx, prepGeom, pt)) {
//                 hash[i] = attempts;
//                 GEOSGeom_destroy_r(ctx, pt);
//                 break;
//             }
//             GEOSGeom_destroy_r(ctx, pt);
//         }
//     }
// 
//     return hash;
// } 


// //selecting filled cells randomly
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds) {
//     vector<int> hash(hashLength, 0);  // Result hash vector
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
//     double cellWidth = (bbox.getXMax() - bbox.getXMin()) / gridSizeX;
//     double cellHeight = (bbox.getYMax() - bbox.getYMin()) / gridSizeY;
// 
//     // Flatten occupied cells into a vector for easy sampling
//     vector<int> occupiedCellList(occupiedCells.begin(), occupiedCells.end());
// 
//     // Step 1: Initialize RNGs
//     vector<std::mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
// 
//     // Step 2: For each hash dimension, sample a random occupied cell
//     for (int i = 0; i < hashLength; ++i) {
//         std::mt19937& generator = generators[i % seeds.size()];
//         std::uniform_int_distribution<int> cellDist(0, static_cast<int>(occupiedCellList.size()) - 1);
// 
//         int cellIndex = occupiedCellList[cellDist(generator)];
//         int gridX = cellIndex % gridSizeX;
//         int gridY = cellIndex / gridSizeX;
// 
//         // Compute dart point as the center of the occupied cell
//         double dartX = bbox.getXMin() + (gridX + 0.5) * cellWidth;
//         double dartY = bbox.getYMin() + (gridY + 0.5) * cellHeight;
// 
//         // Store the number of attempts as 1 (always succeeds)
//         hash[i] = 1;  // or use cellIndex if you want hash to reflect spatial pattern
//     }
// 
//     return hash;
// }

// // Precompute cosine and sine for 360 steps
// const int ANGLE_TABLE_SIZE = 360;
// vector<pair<double, double>> buildTrigTable() {
//     vector<pair<double, double>> table(ANGLE_TABLE_SIZE);
//     for (int i = 0; i < ANGLE_TABLE_SIZE; ++i) {
//         double theta = 2 * M_PI * i / ANGLE_TABLE_SIZE;
//         table[i] = {cos(theta), sin(theta)};
//     }
//     return table;
// }
// //with lookup table for theta
// vector<int> hash2d(const vector<Point>& polygon, int hashLength,
//                    const BoundingBox& bbox,
//                    const vector<vector<GEOSGeometry*>>& grid,
//                    const unordered_set<int>& occupiedCells,
//                    const vector<unsigned int>& seeds) {
//     vector<int> hash(hashLength, 0);
// 
//     int gridSizeX = grid.size();
//     int gridSizeY = grid[0].size();
// 
//     // Step 1: Local polygon MBR
//     double polyXMin = numeric_limits<double>::max(), polyXMax = numeric_limits<double>::lowest();
//     double polyYMin = numeric_limits<double>::max(), polyYMax = numeric_limits<double>::lowest();
// 
//     for (const auto& point : polygon) {
//         polyXMin = min(polyXMin, point.x);
//         polyXMax = max(polyXMax, point.x);
//         polyYMin = min(polyYMin, point.y);
//         polyYMax = max(polyYMax, point.y);
//     }
// 
//     // Step 2: Center of global MBR
//     double centerX = (bbox.getXMin() + bbox.getXMax()) / 2.0;
//     double centerY = (bbox.getYMin() + bbox.getYMax()) / 2.0;
//     double maxRadius = max(bbox.getXMax() - centerX, bbox.getYMax() - centerY);
// 
//     // Step 3: Random Generators
//     vector<mt19937> generators;
//     for (unsigned int seed : seeds) {
//         generators.emplace_back(seed);
//     }
// 
//     uniform_real_distribution<double> uniformDist(0, 1);           // for radius
//     uniform_int_distribution<int> angleIndexDist(0, ANGLE_TABLE_SIZE - 1);  // for angle index
// 
//     // Step 4: Precompute grid cell size inverses
//     double invCellSizeX = gridSizeX / (bbox.getXMax() - bbox.getXMin());
//     double invCellSizeY = gridSizeY / (bbox.getYMax() - bbox.getYMin());
// 
//     // Step 5: Build trig table once
//     static const vector<pair<double, double>> angleTable = buildTrigTable();
// 
//     for (int i = 0; i < hashLength; ++i) {
//         int attempts = 0;
//         mt19937& generator = generators[i % seeds.size()];
// 
//         while (true) {
//             ++attempts;
// 
//             double u = uniformDist(generator);
//             double r = maxRadius * sqrt(u);
// 
//             int angleIdx = angleIndexDist(generator);
//             double cosTheta = angleTable[angleIdx].first;
//             double sinTheta = angleTable[angleIdx].second;
// 
//             double dartX = centerX + r * cosTheta;
//             double dartY = centerY + r * sinTheta;
// 
//             if (dartX < polyXMin || dartX > polyXMax || dartY < polyYMin || dartY > polyYMax)
//                 continue;
// 
//             int gridX = (dartX - bbox.getXMin()) * invCellSizeX;
//             int gridY = (dartY - bbox.getYMin()) * invCellSizeY;
// 
//             if (gridX < 0 || gridX >= gridSizeX || gridY < 0 || gridY >= gridSizeY)
//                 continue;
// 
//             int cellIndex = gridY * gridSizeX + gridX;
// 
//             if (occupiedCells.find(cellIndex) == occupiedCells.end())
//                 continue;
// 
//             if (isPointInsidePolygon({dartX, dartY}, polygon)) {
//                 hash[i] = attempts;
//                 break;
//             }
//         }
//     }
// 
//     return hash;
// }


//generate lsh hash code
vector<int> LSHHashs(vector<double> *sketch, int hashLength, unsigned int seed)
{
    uniform_real_distribution<double> distribution(0.0, static_cast<double>(sketch->size() - 1));
    vector<int> hash = vector<int>(hashLength, 0);

    // Seed the random number generator outside the loop
    minstd_rand generator(seed);

    for (int i = 0; i < hashLength; i++)
    {
        while (true)
        {
            double r = distribution(generator);
            if (r <= int(r) + sketch->at(static_cast<int>(r)))
                break;

            hash[i]++;
        }
       // cout << "Hash[" << i << "] = " << hash[i] << endl;
    }

    return hash;
}



std::vector<std::pair<double, double>> extractCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry)
{
    std::vector<std::pair<double, double>> coordinates;

    int geometryType = GEOSGeomTypeId_r(ctx, geometry);
    if (geometryType != GEOS_POLYGON)
    {
        std::cerr << "Invalid geometry type. Expected Polygon." << std::endl;
        return coordinates;
    }

    // Get the exterior ring of the polygon
    const GEOSGeometry *exteriorRing = GEOSGetExteriorRing_r(ctx, geometry);
    if (exteriorRing == nullptr)
    {
        std::cerr << "Error getting exterior ring." << std::endl;
        return coordinates;
    }

    const GEOSCoordSequence *coordSeq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
    if (coordSeq == nullptr)
    {
        std::cerr << "Error getting coordinate sequence." << std::endl;
        return coordinates;
    }

    unsigned int numPoints;
    GEOSCoordSeq_getSize_r(ctx, coordSeq, &numPoints);

    for (unsigned int i = 0; i < numPoints; ++i)
    {
        double x, y;
        GEOSCoordSeq_getX_r(ctx, coordSeq, i, &x);
        GEOSCoordSeq_getY_r(ctx, coordSeq, i, &y);

        coordinates.push_back(std::make_pair(x, y));
    }

    return coordinates;
    

}


std::pair<std::pair<double, double>, std::pair<double, double>> findBoundingBox(const std::vector<std::pair<double, double>>& points)
{
    if (points.empty()) {
        // Return a default value or handle the case where points are empty
        return {{0.0, 0.0}, {0.0, 0.0}};
    }

    // Initialize the bounding box coordinates
    double xmin = points[0].first;
    double xmax = points[0].first;
    double ymin = points[0].second;
    double ymax = points[0].second;

    // Find the minimum and maximum coordinates
    for (const auto& point : points) {
        xmin = std::min(xmin, point.first);
        xmax = std::max(xmax, point.first);
        ymin = std::min(ymin, point.second);
        ymax = std::max(ymax, point.second);
    }

    // Return the bounding box coordinates
    return {{xmin, ymin}, {xmax, ymax}};
}


GEOSGeometry *createPolygonFromBoundingBox(GEOSContextHandle_t ctx, const std::pair<std::pair<double, double>, std::pair<double, double>> &boundingBox)
{
    // Create a linear ring from the bounding box coordinates
    GEOSCoordSequence *coordSeq = GEOSCoordSeq_create_r(ctx, 5, 2);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 0, boundingBox.first.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 0, boundingBox.first.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 1, boundingBox.second.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 1, boundingBox.first.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 2, boundingBox.second.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 2, boundingBox.second.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 3, boundingBox.first.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 3, boundingBox.second.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 4, boundingBox.first.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 4, boundingBox.first.second);

    // Create a linear ring and add it to a polygon
    GEOSGeometry *linearRing = GEOSGeom_createLinearRing_r(ctx, coordSeq);
    GEOSGeometry *polygon = GEOSGeom_createPolygon_r(ctx, linearRing, nullptr, 0);

    return polygon;
}

//print coordinates of a geometry
void printCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry)
{
    int geometryType = GEOSGeomTypeId_r(ctx, geometry);
    if (geometryType != GEOS_POLYGON)
    {
        cerr << "Invalid geometry type. Expected Polygon." << endl;
        return;
    }

    // Get the exterior ring of the polygon
    const GEOSGeometry *exteriorRing = GEOSGetExteriorRing_r(ctx, geometry);
    if (exteriorRing == nullptr)
    {
        cerr << "Error getting exterior ring." << endl;
        return;
    } 
    const GEOSCoordSequence *coordSeq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
    if (coordSeq == nullptr)
    {
        cerr << "Error getting coordinate sequence." << endl;
        return;
    }
    unsigned int numPoints;
    GEOSCoordSeq_getSize_r(ctx, coordSeq, &numPoints);

    cout << "Number of Points: " << numPoints << endl;
    for (unsigned int i = 0; i < numPoints; ++i)
    {
        double x, y;
        GEOSCoordSeq_getX_r(ctx, coordSeq, i, &x);
        GEOSCoordSeq_getY_r(ctx, coordSeq, i, &y);

        cout << "Point " << i + 1 << ": (" << x << ", " << y << ")" << endl;
    }
   // cout << "error 7" << endl;

    cout << endl;
}