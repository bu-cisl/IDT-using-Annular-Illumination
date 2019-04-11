# High-speed *in vitro* intensity diffraction tomography
MATLAB implementation of paper: **High-speed *in vitro* intensity diffraction tomography**. We provide test data and several quick demo.


### Citation
If you find this project useful in your research, please consider citing our paper:  
[**Jiaji Li, Alex Matlock, Yunzhe Li, Qian Chen, Chao Zuo and Lei Tian, "High-speed *in vitro* intensity diffraction tomography", ariXv preprint](https://arxiv.org)


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
MATLAB R2018a


### How to run the demo
Make sure the data and code at same path and run Main_AI_IDT

[Main_AI_IDT_Diatom_I.m](Main_AI_IDT_Diatom_I.m) for Diatom_I result  
[Main_AI_IDT_Diatom_II.m](Main_AI_IDT_Diatom_II.m) for Diatom_II result  
[Main_AI_IDT_Cheek.m](Main_AI_IDT_Cheek.m) for Cheek cell culster result  
[Main_AI_IDT_Elegans.m](Main_AI_IDT_Elegans.m) for cropped *C. elegans* result  
Functions in [Utilities](Utilities) folder are from [Waller-Lab](https://github.com/Waller-Lab/Angle_SelfCalibration)  


Raw datas are available on Google drive at https://googledirve.com/  


### Illustration of aIDT imaging system
<p align="center">
 <img src="/figs/Imaging system.png">
</p>

### Results  
#### Single cell RI tomography of unstained human cheek cell clusters
<p align="center">
 <img src="/figs/Cheek cell.png" width = "450">
</p>

#### Time-lapse in vitro tomographic imaging of *C. elegans*
<p align="center">
 <img src="/figs/C. elegans.png">
</p>



## License
This project is licensed under the terms of the BSD-3-Clause license. see the [LICENSE](LICENSE) file for details












