# How to use this package #

## Download the package ##

Either clone the repository from github, or download the code. If you downloaded the code as a zip, unzip the file.

## Loading the package ##

All of this first requires [Macaulay2](https://github.com/Macaualy2/M2) installed

### Option 1: Use from the same directory as HHLResolutions.m2 ###

From the same directory as HHLResolutions.m2 (Note: If you want to move the HHLResolutions.m2 file, then you also have to move the HHLResolutions folder with it), run

```
needsPackage "HHLResolutions"
```

### Option 2: Add the directory to the package path ###

If the HHLResolutions.m2 is at the path /home/username/HHLResolutions/HHLResolutions.m2, then run the following:

```
path = append(path, "/home/username/HHLResolutions")
needsPackage "HHLResolutions"
```


### Option 3: Use the FileName optional parameter to needsPackage ###

If the HHLResolutions.m2 is at the path /home/username/HHLResolutions/HHLResolutions.m2, then run the following:

```
needsPackage("HHLResolutions",FileName=>"/home/username/HHLResolutions/HHLResolutions.m2")
```
