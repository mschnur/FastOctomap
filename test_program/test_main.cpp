#include <octomap/octomap.h>
#include <octomap/octomap_timing.h>
#include <string>
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <ios>

#include "FastOctree.h"

#define INSERT_RAY_BY_RAY 0

void printUsage(const char* self)
{
    std::cerr << "\nUSAGE: " << self << " <InputFile>.log/graph\n\n";
    exit(0);
}

void makeOctomapNodeCsv(const std::string& filename, octomap::OcTree& tree)
{
	std::ofstream outFile(filename, std::ios::trunc);
	outFile << "depth,centerX,centerY,centerZ,size,value" << std::endl;
	
	for(octomap::OcTree::tree_iterator it = tree.begin_tree(), end=tree.end_tree(); it!= end; ++it)
	{
		octomap::point3d center = it.getCoordinate();
	    outFile << it.getDepth() << "," << center.x() << "," << center.y() << "," << center.z() << ","
		        << it.getSize() << "," << it->getValue() << std::endl;
	}
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

	std::string extension = logFilename.substr(logFilename.find_last_of(".") + 1);
	std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c){ return std::tolower(c); });
	if(extension == "log") {
		std::cout << "Detected plain ASCII file" << std::endl;
		graph->readPlainASCII(logFilename);
	} else if (extension == "graph") {
		std::cout << "Detected binary graph file" << std::endl;
		graph->readBinary(logFilename);
	} else {
		std::cerr << "Unknown file type" << std::endl;
		return -1;
	}

    std::cout << "\nConstructing OcTree object\n===========================\n";
    std::unique_ptr<octomap::OcTree> tree(new octomap::OcTree(0.05));
    tree->expand();

    std::cout << "\nConstructing FastOctree object\n===========================\n";
    Octree fastOctree = {};
    initOctree(&fastOctree);

    std::cout << "\nAdding Point Clouds\n===========================\n";

    size_t numScans = graph->size();
	size_t currentScan = 1;
	for (octomap::ScanGraph::iterator scan_it = graph->begin(); scan_it != graph->end(); scan_it++) {
		timeval start = {}, stop = {};
		gettimeofday(&start, NULL);  // start timer
#if !INSERT_RAY_BY_RAY
		tree->insertPointCloud((*scan_it)->scan, (*scan_it)->pose.trans());
#endif
        gettimeofday(&stop, NULL);  // start time
#if !INSERT_RAY_BY_RAY
		std::stringstream filenameStream;
		filenameStream << "octomap_nodes_after_pointcloud_" << currentScan << ".csv";
		makeOctomapNodeCsv(filenameStream.str(), *tree);
#endif

		{
			Vector3d sensorOrigin = {};
			const octomap::point3d& sensorOriginOctomap = (*scan_it)->pose.trans();
			initVector3d(&sensorOrigin, sensorOriginOctomap.x(), sensorOriginOctomap.y(), sensorOriginOctomap.z());

			static Vector3d pointsBuffer[300000] = {};
#if INSERT_RAY_BY_RAY
			size_t numPoints = (*scan_it)->scan->size();
#else
			size_t numPoints = 10;
#endif
			for (size_t i = 0; i < numPoints; ++i)
			{
				const octomap::point3d& pt = (*scan_it)->scan->getPoint(i);
				initVector3d(&(pointsBuffer[i]), pt.x(), pt.y(), pt.z());

#if INSERT_RAY_BY_RAY
				tree->insertRay((*scan_it)->pose.trans(), pt);
				std::stringstream filenameStream;
				filenameStream << "octomap_nodes_after_pointcloud_" << currentScan << "_point_" << i << ".csv";
				makeOctomapNodeCsv(filenameStream.str(), *tree);

				std::stringstream fastOctreeFilenameStream;
				fastOctreeFilenameStream << "FastOctree_nodes_after_pointcloud_" << currentScan << "_point_" << i << ".csv";
				insertPointCloud(&fastOctree, &(pointsBuffer[i]), 1, &sensorOrigin);
				createNodeCsv(&fastOctree, fastOctreeFilenameStream.str().c_str());
#endif
			}
#if !INSERT_RAY_BY_RAY
			insertPointCloud(&fastOctree, pointsBuffer, numPoints, &sensorOrigin);

			std::stringstream fastOctreeFilenameStream;
			fastOctreeFilenameStream << "FastOctree_nodes_after_pointcloud_" << currentScan << ".csv";

			createNodeCsv(&fastOctree, fastOctreeFilenameStream.str().c_str());
#endif
		}
		
       
        double time_to_insert = (stop.tv_sec - start.tv_sec) + 1.0e-6 *(stop.tv_usec - start.tv_usec);
	    //std::cout << "Scan #: " << currentScan << std::endl
		//    << "Time to insert: " << time_to_insert << std::endl
		//    << "Memory Usage: " << tree->memoryUsage() << std::endl 
		//	<< "Full memory usage: " << tree->memoryFullGrid() << std::endl;
			
		std::cout << currentScan << "," << (*scan_it)->scan->size() << "," << time_to_insert << std::endl;
       
        currentScan++;
        break;
    }
	
    return 0;
}