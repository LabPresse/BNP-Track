# BNP-Track

This software package includes several Matlab scripts and auxiliary functions, which implement the computational algorithms for the framework BNP-Track described in:
[https://www.nature.com/articles/s41592-024-02349-9](https://www.nature.com/articles/s41592-024-02349-9)

## Prerequisites

Our code requires [Matlab](https://www.mathworks.com/products/matlab.html) and its [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html).

## Downloading BNP-Track

You can download our code by cloning the git repository:

```markdown
git clone https://github.com/LabPresse/BNP-Track.git
```

Alternatively, GitHub allows users to download a repository as a tarball or zipball. Check [here](https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives) for more details.

## Using BNP-Track

1. Once BNP-Track is downloaded (and extracted), open Matlab and navigate to the folder containing all source files.
2. To get a quick start with GUI, please run `BNP_Track_app` in the Matlab command window.
3. If finer control is preferred, please open [`BNP_track_driver.m`](https://github.com/LabPresse/BNP-Track/blob/main/BNP_track_driver.m), directly modify the parameters there, then run this file.
4. Please refer to our manuscript and [this documentation](https://github.com/LabPresse/BNP-Track/blob/main/doc/notes.pdf) for more details.

## Further Assistance

We are working on the implementation and documentation of BNP-Track. You are highly encouraged to contact us at <weiqing1@asu.edu> or <spresse@asu.edu> for help! Issues and pull requests are also very much appreciated!

## Workarounds

**We seek to address issues directly, but this takes time. Therefore, we list some quick workarounds here for these issues.**

1. Currently, BNP-Track only supports reading image stacks that are already stored as a 3D array in Matlab. However, many users have their frames in the TIFF format. To help import the data, we provide [`tiff2mat.m`](https://github.com/LabPresse/BNP-Track/blob/main/tiff2mat.m), which reads an image stack from a TIFF file[^1]. We also offer [`mat2tiff.m`](https://github.com/LabPresse/BNP-Track/blob/main/mat2tiff.m), which is the inverse of `tiff2mat.m`.

[^1]: The orientation of an image may be changed when loaded using `tiff2mat.m`.
