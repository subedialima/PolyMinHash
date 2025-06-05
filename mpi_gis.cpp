#include "parse_geodata.h"
#include "query.h"
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

using namespace std;
vector<string> splitFileMPI(const string filename, MPI_Comm comm) {
    vector<string> result;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Step 1: Open the file and read all lines (only rank 0)
    vector<string> all_lines;
    if (rank == 0) {
        ifstream infile(filename);
        if (!infile.is_open()) {
            cerr << "Error: Unable to open file " << filename << " on rank 0" << endl;
            MPI_Abort(comm, 1);
        }

        string line;
        while (getline(infile, line)) {
            if (!line.empty()) {
                all_lines.push_back(line);
            }
        }
        infile.close();
    }

    // Step 2: Broadcast the number of lines to all ranks
    int total_lines = all_lines.size();
    MPI_Bcast(&total_lines, 1, MPI_INT, 0, comm);

    // Step 3: Scatter the number of lines to read per rank
    int base_lines_per_rank = total_lines / size;
    int remainder = total_lines % size;
    int start_idx = rank * base_lines_per_rank + min(rank, remainder);
    int lines_to_read = base_lines_per_rank + (rank < remainder ? 1 : 0);

    // Step 4: Send the lines to each rank
    vector<string> local_lines(lines_to_read);
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            int send_start = i * base_lines_per_rank + min(i, remainder);
            int send_count = base_lines_per_rank + (i < remainder ? 1 : 0);
            for (int j = 0; j < send_count; ++j) {
                int len = all_lines[send_start + j].size();
                MPI_Send(&len, 1, MPI_INT, i, 0, comm);
                MPI_Send(all_lines[send_start + j].c_str(), len, MPI_CHAR, i, 1, comm);
            }
        }
        // Copy own lines
        for (int i = 0; i < lines_to_read; ++i) {
            local_lines[i] = all_lines[start_idx + i];
        }
    } else {
        for (int i = 0; i < lines_to_read; ++i) {
            int len;
            MPI_Recv(&len, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
            char *buf = new char[len + 1];
            MPI_Recv(buf, len, MPI_CHAR, 0, 1, comm, MPI_STATUS_IGNORE);
            buf[len] = '\0';
            local_lines[i] = string(buf);
            delete[] buf;
        }
    }

    // Step 5: Write local lines to a file
    string temp_dir = "/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/";
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "%stemp_file_%d.txt", temp_dir.c_str(), rank);

    ofstream outfile(temp_filename);
    if (!outfile.is_open()) {
        cerr << "Error: Unable to create temporary file " << temp_filename << " on rank " << rank << endl;
        MPI_Abort(comm, 1);
    }

    for (const string &line : local_lines) {
        outfile << line << "\n";
    }
    outfile.close();

    result.push_back(temp_filename);
    return result;
}
