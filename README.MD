# High-speed *in vitro* intensity diffraction tomography
MATLAB implementation of paper: **High-speed in vitro intensity diffraction tomography**. We provide test data and a quick demo.


### Citation
If you find this project useful in your research, please consider citing our paper:

[**Jiaji Li, Alex Matlock, Yunzhe Li, Qian Chen, Chao Zuo and Lei Tian, "High-speed in vitro intensity diffraction tomography", ariXv preprint ]


### Abstract
We demonstrate a label-free, scan-free intensity diffraction tomography technique utilizing annular illumination (aIDT) to rapidly characterize large-volume 3D refractive index distributions \emph{in vitro}.
By optimally matching the illumination geometry to the microscope pupil, our technique reduces the data requirement by 60$\times$ to achieve high-speed 10 Hz volume rates.
Using 8 intensity images, we recover \mbox{$\sim350\times100\times20\mu$m$^3$} volumes with near diffraction-limited lateral resolution of 487 nm and axial resolution of 3.4 $\mu$m.
Our technique's large volume rate and high resolution enables 3D quantitative phase imaging of complex living biological samples across multiple length scales.
We demonstrate aIDT's capabilities on unicellular diatom microalgae, epithelial buccal cell clusters with native bacteria, and live \emph{Caenorhabditis elegans} (\emph{C. elegans}) specimens.
Within these samples, we recover macroscale cellular structures, subcellular organelles, and dynamic micro-organism tissues with minimal motion artifacts.
Quantifying such features has significant utility in oncology, immunology, and cellular pathophysiology, where these morphological features are evaluated for changes in the presence of disease, parasites, and new drug treatments.
aIDT shows promise as a powerful high-speed, label-free microscopy technique for these applications where natural imaging is required to evaluate environmental effects on a sample in real-time.


### Requirements
MATLAB R2018a


### Illustration of our aIDT imaging system
<p align="center">
 <img src="/figs/Imaging system.png">
</p>


### How to run the demo
Make sure you have all dependent packages installed correctly and run [demo.py](demo.py).
To reproduce reconstructions run the corresponding MATLAB scripts as described below.

Main_AI_IDT_Diatom_I.m
Main_AI_IDT_Diatom_II.m
Main_AI_IDT_Cheek.m

Raw datas are available on Google drive at https://googledirve.com/
3 Frame Elegans data contain raw captured images

### Results
## Single cell RI tomography of unstained human cheek cell clusters
<p align="center">
 <img src="/figs/Cheek cell.png">
</p>

## On UNSEEN cell type
<p align="center">
 <img src="/figs/C. elegans.png">
</p>



## License
This project is licensed under the terms of the BSD-3-Clause license. see the [LICENSE](LICENSE) file for details












