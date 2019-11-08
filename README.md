# FastOctomap
Just once, you'll need to extract the large data txt file, which can be done by running the script
```
./prepare_data_file.sh
```

Build with
```
./build.sh
```

Currently the octomap sub-repository has a branch mps/revelles_ray_trace. This branch is currently a work in progress,
so the octomap sub-repository might fail when running build.sh.

Once built, the test programs can be run with
```
./deploy/original_octomap_test/bin/test_original_octomap rgbd_dataset_freiburg1_desk2.txt
./deploy/octomap_test/bin/test_mod_octomap rgbd_dataset_freiburg1_desk2.txt
```

The second of those two commands won't work if the octomap sub-repository can't build successfully.