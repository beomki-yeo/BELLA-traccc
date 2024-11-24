## Repository for BELLA experiment reconstruction with traccc

### Prerequisites
- gcc >= 11.2 (C++ 20 support)
- CMake >= 3.30.2
- Boost >= 1.86.0 (program_options, filesystem, log)

### How to run

After compiling the repository, go to the `shell` directory and run the script.
Make sure that `${BUILD_DIR}` in the script is consistent with the build directory

```sh
cd shell
bash telescope_script.sh
```

Then it will create some output csv files at `data` directory
