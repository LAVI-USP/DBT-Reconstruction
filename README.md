[![Python Version](https://upload.wikimedia.org/wikipedia/commons/f/fc/Blue_Python_3.7_Shield_Badge.svg)](https://github.com/LAVI-USP/pyDBT) [![View LAVI DBT-Reconstruction on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/71138-lavi-dbt-reconstruction)   [![Open in Code Ocean](https://codeocean.com/codeocean-assets/badge/open-in-code-ocean.svg)](https://codeocean.com/capsule/1859673/) 

Open source reconstruction toolbox for digital breast tomosynthesis
======

**Objective:** Open-source reconstruction toolbox, which features DBT reconstruction with adjustable acquisition geometries. This toolbox is intended for academic usage and it aims at expanding the academic research in this field. Since it is an open-source toolbox, researchers are welcome to contribute to the current version of the software.

## Abstract:

 Digital breast tomosynthesis (DBT) is an emerging imaging modality used for breast cancer screening. This method minimizes the limitation of tissue superimposition from conventional digital mammography. The reconstruction of the 3D volume from radiographic projections is an important step and it is still a challenging task. Clinical units from various manufacturers differ in acquisition geometry, and thus different reconstruction parameters must be set for each of them. We have developed an open-source MATLAB toolbox for DBT reconstruction. The proposed toolbox is intended for academic usage. The virtual Shepp-Logan phantom and a physical BR3D phantom were used to validate the proposed implementation of the reconstruction methods and the acquisition geometry of the clinical unit. Visual assessment of the reconstructed slices obtained with the proposed toolbox indicates that the reconstruction was performed successfully. The results were evaluated using a structural similarity metric. Visual comparison against the slices obtained with the built-in commercial reconstruction software shows the potential of the proposed reconstruction pipeline, which can be easily modified to incorporate sophisticated algorithms and accelerated using parallel computation techniques.

## How to install?

 1. Download the toolbox or clone the directory [here](https://github.com/LAVI-USP/DBT-Reconstruction.git).
 
 2. The toolbox was tested on **Windows** 10 x64, **Linux** (Ubuntu) x64, and **MacoS** High Sierra machines, with Matlab 2015a and 2017b.

 3. In order to use CUDA or mex files, please follow the instructions [here](https://github.com/LAVI-USP/DBT-Reconstruction/blob/master/Functions/Sources/README.md). These files were only compiled and tested on **Windows** 10 x64.
 
 4. If you don't want to bother with these steps, you might want to run the toolbox on [Codeocean](https://codeocean.com/capsule/1859673/).

 5. Finally, you can find the DBT toolbox **Python version** [here](https://github.com/LAVI-USP/pyDBT). 

Please report issues [here](https://github.com/LAVI-USP/DBT-Reconstruction/issues).

## Contribute? 

We are pleased with any contributions. Feel free to make any [pull requests](https://github.com/LAVI-USP/DBT-Reconstruction/pulls) or send us an e-mail.


## Some results:

**DBT projection simulation:**

![Imgur](https://i.imgur.com/tcgqs8I.gif?1)


**Projections:**

![Imgur](https://i.imgur.com/VXk6MFi.gif?2)


**Reconstructed volume:**

![Imgur](https://i.imgur.com/R9tB35f.gif?2)

**Speedup with CUDA and OpenMP implementations:**

![Imgur](https://i.imgur.com/PreHXUz.png)

It was reconstructed a volume of 1978x1058x107 voxels. The DD stands for Distance-driven and DDb for Branchless-Distance-driven. 

## Toolbox manual:

I created a [wiki page](https://github.com/LAVI-USP/DBT-Reconstruction/wiki/Toolbox-Manual) with the manual. Please, send me any advice or correction related to this file.

## Contact:

If you have any questions or suggestions, please send us an e-mail:

- Rodrigo - rodrigo dot vimieiro at gmail dot com
- Marcelo - mvieira at sc dot usp dot br

## License:

The toolbox is licensed under the **GNU General Public License v3.0**. Please check the [licence file](https://github.com/LAVI-USP/DBT-Reconstruction/blob/master/LICENSE).

## Reference:

If you use the toolbox, we will be very grateful if you refer to this [paper](https://doi.org/10.1007/978-981-13-2517-5_53):

> Vimieiro R.B., Borges L.R., Vieira M.A.C. (2019) Open-Source Reconstruction Toolbox for Digital Breast Tomosynthesis. In: Costa-Felix R., Machado J., Alvarenga A. (eds) XXVI Brazilian Congress on Biomedical Engineering. IFMBE Proceedings, vol 70/2. Springer, Singapore.

## Citations:

You can find [here](https://scholar.google.com.br/scholar?oi=bibs&hl=pt-BR&cites=3156269064066227282) the papers that have used the toolbox.

## Acknowledgments:

This work was supported by the São Paulo Research Foundation ([FAPESP](http://www.fapesp.br/) grant 2016/25750-0) and by the National Council for Scientific and Technological Development ([CNPq](http://www.cnpq.br/)). Nobody does anything alone, so we would like to thank the contribution of our lab members and the [Barretos Love Hospital](https://www.hcancerbarretos.com.br) for providing the images of DBT.

---

Laboratory of Computer Vision ([Lavi](http://iris.sel.eesc.usp.br/lavi/))  
Department of Electrical and Computer Engineering  
São Carlos School of Engineering, University of São Paulo  
São Carlos - Brazil
