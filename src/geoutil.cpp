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
                
                
            }
        }
    }

    result.shrink_to_fit();
    return result;
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




// based on global MBR 
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
