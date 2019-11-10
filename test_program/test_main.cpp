#include <octomap/octomap.h>
#include <octomap/octomap_timing.h>
#include <string>
#include <iostream>
#include <memory>

void printUsage(const char* self)
{
    std::cerr << "\nUSAGE: " << self << " <InputFile>.log\n\n";
    exit(0);
}

int main(int argc, char** argv)
{
    // default values:
    std::string logFilename = "";

    if (argc != 2) {
        printUsage(argv[0]);
    } else {
        logFilename = std::string(argv[1]);
    }

    std::cout << "\nReading Log file\n===========================\n";
    std::unique_ptr<octomap::ScanGraph> graph(new octomap::ScanGraph());
    graph->readPlainASCII(logFilename);
   
    std::cout << "\nConstructing OcTree object\n===========================\n";
    std::unique_ptr<octomap::OcTree> tree(new octomap::OcTree(0.05));
    tree->expand();
    

    std::cout << "\nAdding Point Clouds\n===========================\n";

    size_t numScans = graph->size();
	size_t currentScan = 1;
	for (octomap::ScanGraph::iterator scan_it = graph->begin(); scan_it != graph->end(); scan_it++) {
		timeval start = {}, stop = {};
		gettimeofday(&start, NULL);  // start timer
		tree->insertPointCloud((*scan_it)->scan, (*scan_it)->pose.trans());
        gettimeofday(&stop, NULL);  // start time
       
        double time_to_insert = (stop.tv_sec - start.tv_sec) + 1.0e-6 *(stop.tv_usec - start.tv_usec);
	    //std::cout << "Scan #: " << currentScan << std::endl
		//    << "Time to insert: " << time_to_insert << std::endl
		//    << "Memory Usage: " << tree->memoryUsage() << std::endl 
		//	<< "Full memory usage: " << tree->memoryFullGrid() << std::endl;
			
		std::cout << currentScan << "," << (*scan_it)->scan->size() << "," << time_to_insert << std::endl;
       
        currentScan++;
    }
	
    return 0;
}