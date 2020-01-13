overlapping community detection on func conn mats
==============================

Check out the scripts in this repo [here](https://github.com/faskowit/overlap_on_NxN_func/tree/master/src/scripts) for running some overlapping community detection algorithms on thresholded func mats. You'll have to install the actual algorithms appropriately on your computer. Algos:
* [SVINET](https://github.com/premgopalan/svinet)
* [AGMfit](http://snap.stanford.edu/agm/)
* [commDetNMF](https://github.com/ipsorakis/commDetNMF)
* and we made our own, link clustering, based on weighted line graphs (ThrLink)

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── bin
    ├── config
    ├── data
    │   ├── external
    │   ├── interim
    │   ├── processed
    │   └── raw
    ├── docs
    ├── LICENSE
    ├── README.md
    ├── reports
    │   └── figures
    └── src
        ├── analysis
        ├── external
        ├── funcs
        ├── processing
        ├── scripts
        └── viz
