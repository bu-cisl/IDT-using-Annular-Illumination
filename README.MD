# High-speed *in vitro* intensity diffraction tomography using Annular Illumination (aIDT)
MATLAB implementation of paper: **High-speed *in vitro* intensity diffraction tomography**. This github provides the codes, links to test data, and several quick demos illustrating this technique.


### Citation
If you use this project in your own research, please cite our corresponding paper: 
[Jiaji Li, Alex Matlock, Yunzhe Li, Qian Chen, Chao Zuo and Lei Tian, "High-speed *in vitro* intensity diffraction tomography", arXiv preprint](https://arxiv.org/abs/1904.06004)


### Abstract
We demonstrate a label-free, scan-free intensity diffraction tomography technique utilizing annular illumination (aIDT) to rapidly characterize large-volume 3D refractive index distributions *in vitro*.
By optimally matching the illumination geometry to the microscope pupil, our technique reduces the data requirement by 60X to achieve high-speed 10 Hz volume rates.
Using 8 intensity images, we recover ~ 350um X 100um X 20um volumes with near diffraction-limited lateral resolution of 487 nm and axial resolution of 3.4 um.
Our technique's large volume rate and high resolution enables 3D quantitative phase imaging of complex living biological samples across multiple length scales.
We demonstrate aIDT's capabilities on unicellular diatom microalgae, epithelial buccal cell clusters with native bacteria, and live *Caenorhabditis elegans* (*C. elegans*) specimens.
Within these samples, we recover macroscale cellular structures, subcellular organelles, and dynamic micro-organism tissues with minimal motion artifacts.
Quantifying such features has significant utility in oncology, immunology, and cellular pathophysiology, where these morphological features are evaluated for changes in the presence of disease, parasites, and new drug treatments.
aIDT shows promise as a powerful high-speed, label-free microscopy technique for these applications where natural imaging is required to evaluate environmental effects on a sample in real-time.


### Requirements
This code was developed using MATLAB R2018a. Please be advised that this code may not work with other editions of MATLAB.

### How to run the demo
1. Download this Repository to your computer.
2. Download the raw data from our [Google drive](https://drive.google.com/drive/folders/1E8U1YWQ14bDbgwMJnaUZl1I6Lq3XoGWj?usp=sharing). 
3. Make sure the data and code are in the same folder path and run a "Main_AI_IDT" script. The following script versions are available for each demo:  
-[Main_AI_IDT_Diatom_I.m](Main_AI_IDT_Diatom_I.m) for Diatom_I result  
-[Main_AI_IDT_Diatom_II.m](Main_AI_IDT_Diatom_II.m) for Diatom_II result  
-[Main_AI_IDT_Cheek.m](Main_AI_IDT_Cheek.m) for Cheek cell culster result  
-[Main_AI_IDT_Elegans.m](Main_AI_IDT_Elegans.m) for cropped *C. elegans* result  

We have implemented the functions in [Utilities](Utilities) from [Professor Laura Waller's Lab](https://github.com/Waller-Lab/Angle_SelfCalibration)  

### Illustration of aIDT imaging system
<p align="center">
 <img src="/figs/Imaging system.png">
</p>

### Results  
#### Single cell RI tomography of unstained human cheek cell clusters
<p align="center">
 <img src="/figs/Cheek cell.png" width = "450">
</p>

#### Time-lapse *in vitro* tomographic imaging of *C. elegans*
<p align="center">
 <img src="/figs/C. elegans.png">
</p>



## License
This project is licensed under the terms of the BSD-3-Clause license. see the [LICENSE](LICENSE) file for details












