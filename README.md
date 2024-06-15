# GPSCO

## GPSCO: Global Planar Structural Constraint Optimal based point cloud registration algorithm for repetitive structures

![flowchart.jpg](flowchart.jpg)

---

### Build

#### Requirement

 ```
    cmake(>=3.16)
    PCL(>=1.10)
    Ceres(2.1.0)
    spdlog(1.12.0)
 ```

#### Option 1(shell)

 ```
    cd ${GPSCO_DIR}
    mkdir build
    cd build
    cmake ..
    make -j14
 ```
*Note:* If the number of threads on your computer is less than `16`, 
modify `make -j14` to match the number of threads on your computer, such as `make -j4`.

#### Option 2(IDE)

Use any IDE that can directly handle CMakeLists files to open the CMakeLists.txt in the root directory of GPSCO. I recommend using [CLion](https://www.jetbrains.com/clion/).

---

### Usage

* register a 'source' point cloud to a 'target' point cloud.

```
    ./examples/exam_regis source.ply target.ply result.txt
```

---

### Benchmark Dataset

HS([Download Link1](https://drive.google.com/drive/folders/1OFHm4iSt0wIle2MeTROb93jYmVoc37mR?usp=sharing) And [Download Link2](https://pan.baidu.com/s/14ZvvR6qfYkQkbzvjEzspjg?pwd=lwcl)): Real-world Scans with Repetitive Structures. This dataset is part of the GPSCO work.

*Note:* The HS dataset was created to complement the open source dataset in terms of repetitive structure. 
The dataset includes raw scans, ground truth, and point cloud overlap.

### Acknowledgments

We thank the respective authors for open sourcing their methods.

* [SC2-PCR](https://github.com/ZhiChen902/SC2-PCR)
* [Super4PCS](https://github.com/nmellado/Super4PCS)
* [PLADE](https://github.com/chsl/PLADE)
* [VPFBR](https://github.com/zhanjiawang/VPFBR-L)
