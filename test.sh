#!/bin/bash

rm *.csv
rm -rf ./build/ ./deploy/

./build.sh fast

./deploy/octomap_test/bin/test_mod_octomap rgbd_dataset_freiburg1_desk2.graph


prefix=""
if [ "$1" == "conv" ]; then
    prefix="conv_"
fi

echo "Point 1"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_0.csv "$prefix"octomap_nodes_after_pointcloud_1_point_0.csv

echo "Point 2"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_1.csv "$prefix"octomap_nodes_after_pointcloud_1_point_1.csv

echo "Point 3"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_2.csv "$prefix"octomap_nodes_after_pointcloud_1_point_2.csv

echo "Point 4"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_3.csv "$prefix"octomap_nodes_after_pointcloud_1_point_3.csv

echo "Point 5"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_4.csv "$prefix"octomap_nodes_after_pointcloud_1_point_4.csv

echo "Point 6"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_5.csv "$prefix"octomap_nodes_after_pointcloud_1_point_5.csv

echo "Point 7"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_6.csv "$prefix"octomap_nodes_after_pointcloud_1_point_6.csv

echo "Point 8"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_7.csv "$prefix"octomap_nodes_after_pointcloud_1_point_7.csv

echo "Point 9"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_8.csv "$prefix"octomap_nodes_after_pointcloud_1_point_8.csv

echo "Point 10"
diff "$prefix"FastOctree_nodes_after_pointcloud_1_point_9.csv "$prefix"octomap_nodes_after_pointcloud_1_point_9.csv