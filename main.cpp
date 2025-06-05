
#define GEOS_USE_ONLY_R_API
#include "parse_geodata.h"
#include "query.h"
#include "thread.h"
#include <cmath>
#include <list>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map> 
#include <set>
#include <geos_c.h>
#include "mpi_gis.h" 
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//typedef std::unordered_map<std::string, QueryResult> AllResults;

using namespace std;

//Define some data types for readability
using ResultType = std::vector<std::pair<std::string, std::vector<std::pair<double, std::string>>>>;

using NeighborPair = std::pair<double, std::string>; // Represents distance and Neighbor ID

using QueryResult = std::vector<NeighborPair>; // All neighbors for a query

// Helper function to check if a file exists
bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}




//prasing output from a file 
void ParseFromFile(GEOSContextHandle_t ctx, const std::string& filename, std::unordered_map<std::string, QueryResult>& allResults) {
    std::ifstream file(filename);
    std::string line;
    std::string currentQueryID;

    while (std::getline(file, line)) {
        // Check for Query ID line
        if (line.find("Query ID:") != std::string::npos) {
            // Extract the Query ID
            std::string idPrefix = "Query ID: ";
            size_t idStartPos = line.find(idPrefix) + idPrefix.length();
            currentQueryID = line.substr(idStartPos);
        }
        // Check for Neighbor line
        else if (!currentQueryID.empty() && line.find("Neighbor ID:") != std::string::npos) {
            std::string neighborPrefix = "Neighbor ID: ";
            size_t neighborStartPos = line.find(neighborPrefix) + neighborPrefix.length();
            size_t commaPos = line.find(", Distance:");
            std::string neighborID = line.substr(neighborStartPos, commaPos - neighborStartPos);

            std::string distancePrefix = "Distance: ";
            size_t distanceStartPos = line.find(distancePrefix) + distancePrefix.length();
            double distance = std::stod(line.substr(distanceStartPos));

            // Add the neighbor information to the current query's list
            allResults[currentQueryID].emplace_back(distance, neighborID);
        }
        // Reset currentQueryID if the line is empty, indicating a new query section
        else if (line.empty()) {
            currentQueryID.clear();
        }
    }

    file.close();
}




int main(int argc, char **argv)
{
    int rank, size;
    int numSplit;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm = MPI_COMM_WORLD;
    char hostname[256];
    gethostname(hostname, 256);

    if (argc != 4)
    {
        if (rank == 0)
        {
            cerr << "Incorrect number of arguments!\n";
            cerr << "Usage: mpirun -n <num_processes> ./parallel_main <dataFile> <inputFile>\n";
        }
        MPI_Finalize();
        return 1;
    }
    string datafile = string(argv[1]);
    string inputfile = string(argv[2]);
    unsigned int numNeighbor = atoi(argv[3]);
    
    


    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    GEOSWKTWriter* wktWriter = GEOSWKTWriter_create_r(ctx);

    vector<string> localFiles = splitFileMPI(datafile, comm);
    
   // vector<GEOSGeometry *> *geos = read_wkt(localFiles[0], ctx); 
    unordered_map<string, GEOSGeometry*> *geos = read_wkt_with_id(localFiles[0], ctx); //need to store unique id as well
    
    
    //print id + geometries
    if (geos != nullptr && !geos->empty()){
        for (const auto& pair : *geos) {
            const string& id = pair.first;
            GEOSGeometry* geometry = pair.second;

            // cout << "ID: " << id << " - input Geometry: ";
            // printCoordinates(ctx, geometry);
            cout << endl;
        }
    }
    else {
        cerr << "No geometries to print or failed to read geometries." << endl;
    }

    if (geos->empty())
    {
        cerr << "Error: No geometries to process in process " << rank << endl;
        MPI_Finalize();
        GEOS_finish_r(ctx);
        return 1;
    }
    
    //read query geometries and their unique id
    // vector<string> localFiles = splitFileMPI(datafile, comm);
    unordered_map<string, GEOSGeometry*> *querygeos = read_wkt_with_id(inputfile, ctx); 

    // Vector to hold centered geometries
    unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID;

    // Iterate over the unordered_map, center each geometry, and store it with the identifier
    for (const auto& pair : *querygeos) {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            // Store the centered geometry with its identifier in the map
            centeredQueryGeosWithID[pair.first] = centeredGeo;
        } else {
            std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
        }
    }

    //perform lsh query
    int gridSize = 25;
    int hashLength = 1;
    ResultType LSHoutput;

    //  LSHHashFunction(*geos, rank, centeredqueryGeos, numNeighbor, gridSize, hashLength, LSHoutput, ctx);
  
    //calling lsh function which performs lsh query process   
    LSHHashFunction(*geos, rank, centeredQueryGeosWithID, numNeighbor, gridSize, hashLength, LSHoutput, ctx);

    //for writing output in a file
    std::string outputFilename = "LSH_results_rank_" + std::to_string(rank) + ".txt";
    std::ofstream outFile("/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/results/LSH_results_rank_" + std::to_string(rank) + ".txt", std::ios::out);

    if (!outFile.is_open()) {
        std::cerr << "Failed to open file for writing results." << std::endl;
        MPI_Finalize();
        return 1; // or handle the error appropriately
    }

    // Assuming 'outFile' is successfully opened as per your code

    for (const auto& queryResult : LSHoutput) {
        const std::string& queryID = queryResult.first;  // Query ID
        const auto& neighbors = queryResult.second;  // Vector of (distance, neighbor ID) pairs

        // Write the query ID
        outFile << "Query ID: " << queryID << std::endl;

        // Write each neighbor ID and distance
        for (const auto& neighbor : neighbors) {
            outFile << "Neighbor ID: " << neighbor.second << ", Distance: " << neighbor.first << std::endl;
        }

        outFile << std::endl; // Blank line for separating entries
    }

    outFile.close(); // Don't forget to close the file after writing

    //Barrier to synchronize all processes before reading files
    MPI_Barrier(MPI_COMM_WORLD);

    std::unordered_map<std::string, QueryResult> globalResults;
    if (rank == 0) {
        std::string baseDirectory = "/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/results/";
        std::string baseFilename = "LSH_results_rank_";
        std::string fileExtension = ".txt";

       std::unordered_map<std::string, QueryResult> allResult;

        // Process each file individually
        for (int i = 0; i < size; ++i) {
            std::string filename = baseDirectory + baseFilename + std::to_string(i) + fileExtension;
            std::ifstream inFile(filename);
            
            if (inFile.is_open()) {
                //  std::cout << "Reading contents from file: " << filename << std::endl;
                
                // Reset allResult for each file
                allResult.clear();
                
                // ParsefromFile function fills allResult with data from the current file
                ParseFromFile(ctx, filename, allResult);
                
                // Merge fileResults into globalResults
                for (const auto& [queryID, neighbors] : allResult) {
                    // Insert new neighbors and existing ones into a temporary vector
                    std::vector<NeighborPair> combinedNeighbors = globalResults[queryID];
                    combinedNeighbors.insert(combinedNeighbors.end(), neighbors.begin(), neighbors.end());

                    // Sort combined neighbors by distance
                    std::sort(combinedNeighbors.begin(), combinedNeighbors.end(),
                          [](const NeighborPair& a, const NeighborPair& b) { return a.first < b.first; });

                    // Keep only the top k neighbors
                    if (combinedNeighbors.size() > numNeighbor) {
                        combinedNeighbors.resize(numNeighbor);
                    }

                    // Update globalResults with the top 3 neighbors
                    globalResults[queryID] = combinedNeighbors;
                }

                inFile.close();

            } else {
                std::cerr << "Failed to open file: " << filename << std::endl;
            }           
        }
    }



    if (rank == 0) {
        // --- Compute “how close” LSH @1 really is ---
        double sumS1_LSH = 0.0;
        std::vector<double> allS1_LSH;

        for (auto const& [qID, lshNeighbors] : globalResults) {
            if (lshNeighbors.empty()) continue;
            // take the very top-1 neighbor
            double lshDist1 = lshNeighbors[0].first;
            double s1 = 1.0 - lshDist1;               // convert dist → Jaccard similarity
            sumS1_LSH += s1;
            allS1_LSH.push_back(s1);
        }

        int N = allS1_LSH.size();
        if (N > 0) {
            double meanS1 = sumS1_LSH / N;
            std::sort(allS1_LSH.begin(), allS1_LSH.end());
            double medianS1 = allS1_LSH[N/2];
            double p10     = allS1_LSH[N/10];
            double p90     = allS1_LSH[(9*N)/10];
            int cnt80     = std::count_if(allS1_LSH.begin(), allS1_LSH.end(), [](double s){ return s >= 0.8; });

            std::cout << "\n=== LSH Recall@1 True‐Neighbor Similarity Stats ===\n";
            std::cout << "  Queries evaluated: " << N << "\n";
            std::cout << "  Mean similarity     = " << meanS1 << "\n";
            std::cout << "  Median similarity   = " << medianS1 << "\n";
            std::cout << "  10–90 pctile range = [" << p10 << ", " << p90 << "]\n";
            std::cout << "  # with s>=0.8       = " << cnt80 << " ("  << (100.0 * cnt80 / N) << "%)\n\n";
            } else {
                std::cout << "No LSH results found to compute Recall@1 similarity.\n";
            }
    }

    //start brute fore search

    vector<pair<string, vector<pair<double, string>>>> queryBFResults;
    BFquery(rank, centeredQueryGeosWithID, *geos, numNeighbor, queryBFResults) ;
 
    // Construct the output filename using the process rank
    std::string outputBFFilename = "BFquery_results_rank_" + std::to_string(rank) + ".txt";
    // Open an output file stream
    std::ofstream outBFFile("/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/results/BFquery_results_rank_" + std::to_string(rank) + ".txt", std::ios::out);
 
    if (!outBFFile.is_open())
    {
        std::cerr << "Failed to open file for writing results." << std::endl;
        MPI_Finalize();
        return 1; // or handle the error appropriately
    }
     
    for (const auto& resultPair : queryBFResults) {
        const string& queryID = resultPair.first;
        const auto& neighborPairs = resultPair.second;

        // Write query ID instead of geometry WKT
        outBFFile << "Query ID: " << queryID << std::endl;

        for (const auto& neighborPair : neighborPairs) {
            double distance = neighborPair.first;
            const string& neighborID = neighborPair.second;

            // Write neighbor ID instead of geometry WKT
            outBFFile << "Neighbor ID: " << neighborID << ", Distance: " << distance << std::endl;
        }
        outBFFile << std::endl;// Separator for readability; // Separator for readability
    }
    
    
    outBFFile.close(); 

    std::unordered_map<std::string, QueryResult> globalBFResults;
    if (rank == 0) {
        std::string baseDirectory = "/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/results/";
        std::string baseFilename = "BFquery_results_rank_";
        std::string fileExtension = ".txt";

       std::unordered_map<std::string, QueryResult> allBFResult;


        // Process each file individually
        for (int i = 0; i < size; ++i) {
            std::string filename = baseDirectory + baseFilename + std::to_string(i) + fileExtension;
            std::ifstream inFile(filename);
            
            if (inFile.is_open()) {
               // std::cout << "Reading contents from file: " << filename << std::endl;
                
                // Reset allResult for each file
                allBFResult.clear();
                
                // ParsefromFile function fills allResult with data from the current file
                ParseFromFile(ctx, filename, allBFResult);
                
                // Merge fileResults into globalResults
                // Merge fileResults into globalResults
                for (const auto& [queryID, neighbors] : allBFResult) {
                    // Insert new neighbors and existing ones into a temporary vector
                    std::vector<NeighborPair> combinedNeighbors = globalBFResults[queryID];
                    combinedNeighbors.insert(combinedNeighbors.end(), neighbors.begin(), neighbors.end());

                    // Sort combined neighbors by distance
                    std::sort(combinedNeighbors.begin(), combinedNeighbors.end(),
                          [](const NeighborPair& a, const NeighborPair& b) { return a.first < b.first; });

                    // Keep only the top 3 neighbors
                    if (combinedNeighbors.size() > numNeighbor) {
                        combinedNeighbors.resize(numNeighbor);
                    }

                    // Update globalResults with the top 3 neighbors
                    globalBFResults[queryID] = combinedNeighbors;
                }

                inFile.close();

            } 
        }

    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // start comparing performance (currently used)
    if (rank == 0) {

        FILE* output = fopen("/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/comparisonResults.csv", "w");

        if (output == NULL) {
            cerr << "Failed to open file for writing.\n";
        } else {
            fprintf(output, "query_id,avg_recall,mse,false_positive_rate,total_lsh_neighbors,total_bf_neighbors,total_matched_neighbors\n");

            double totalOverallRecall = 0.0;
            double totalOverallBFRecall = 0.0;
            double totalMSE = 0.0;
            int totalCount = 0;
            int totalFalsePositives = 0;
            int totalBFNeighborCount = 0;  // Counter for all BF neighbors, including duplicates
            int totalLSHNeighborCount = 0; // Counter for all LSH neighbors, including duplicates
            int totalMatchedNeighborsCount = 0;

            for (const auto& [queryID, bfNeighbors] : globalBFResults) {
                auto lshResultIt = globalResults.find(queryID);
                if (lshResultIt == globalResults.end()) continue;
                auto& lshNeighbors = lshResultIt->second;
                int correctNeighbors = 0;
                double mse = 0.0;
                int mseidx = 0;

                // Count all BF neighbors (including duplicates)
                totalBFNeighborCount += bfNeighbors.size(); // Add the number of neighbors for this query

                // Count all LSH neighbors (including duplicates)
                totalLSHNeighborCount += lshNeighbors.size(); // Add the number of neighbors for this query

                // Matching BF and LSH neighbors for comparison
                for (const auto& [distance, lshNeighborID] : lshNeighbors) {
                    bool isMatched = false;
                    for (const auto& [bfDistance, bfNeighborID] : bfNeighbors) {
                        if (lshNeighborID == bfNeighborID) {
                            // This neighbor is correct
                            correctNeighbors++;
                            isMatched = true;

                            // Calculate MSE for matching neighbors
                            mse += pow(bfDistance - distance, 2);
                            mseidx++;
                            break;
                        }
                    }

                    // False positive if no match is found in BF results
                    if (!isMatched) {
                        totalFalsePositives++;
                    }
                }

                if (!lshNeighbors.empty()) {
                    double queryRecall = static_cast<double>(correctNeighbors) / lshNeighbors.size();
                    totalOverallRecall += queryRecall;
                }

                if (!bfNeighbors.empty()) {
                    double queryBFRecall = static_cast<double>(correctNeighbors) / bfNeighbors.size();
                    totalOverallBFRecall += queryBFRecall;
                }

                if (mseidx > 0) {
                    mse /= mseidx;
                    totalMSE += mse;
                }

                totalCount++;
                totalMatchedNeighborsCount += correctNeighbors;
            }

            // Calculate averages
            double averageRecall = totalCount > 0 ? totalOverallRecall / totalCount : 0;
            double averageBFRecall = totalCount > 0 ? totalOverallBFRecall / totalCount : 0;
            double averageMSE = totalCount > 0 ? totalMSE / totalCount : 0;
            double falsePositiveRate = totalLSHNeighborCount > 0 ? static_cast<double>(totalFalsePositives) / totalLSHNeighborCount : 0;

            // Output the results
            cout << "Average Recall: " << averageRecall << endl;
            cout << "Average Recall by BF neighbors: " << averageBFRecall << endl;
            cout << "MSE: " << averageMSE << endl;
            cout << "Total Matched Neighbors: " << totalMatchedNeighborsCount << endl;
            cout << "Total LSH Neighbors (including duplicates): " << totalLSHNeighborCount << endl;
            cout << "Total BF Neighbors (including duplicates): " << totalBFNeighborCount << endl;
            cout << "Total False Positives rate: " << falsePositiveRate << endl;

            // Additional console output...

            fclose(output);
        }
    }

    MPI_Finalize();

    return 0;
    
    GEOS_finish_r(ctx);
    MPI_Finalize();
}
