Open source reconstruction toolbox for digital breast tomosynthesis
======

**Objective:** Open-source reconstruction toolbox, which features DBT reconstruction with adjustable acquisition geometries. This toolbox is intended for academic usage and it aims at expanding the academic research in this field. Since it is an open-source toolbox, researchers are welcome to contribute to the current version of the software.

## Abstract:

 Digital breast tomosynthesis (DBT) is an emerging imaging modality used for breast cancer screening. This method minimizes the limitation of tissue superimposition from conventional digital mammography. The reconstruction of the 3D volume from radiographic projections is an important step and it is still a challenging task. Clinical units from various manufacturers differ in acquisition geometry, and thus different reconstruction parameters must be set for each of them. We have developed an open-source MATLAB toolbox for DBT reconstruction. The proposed toolbox is intended for academic usage. The virtual Shepp-Logan phantom and a physical BR3D phantom were used to validate the proposed implementation of the reconstruction methods and the acquisition geometry of the clinical unit. Visual assessment of the reconstructed slices obtained with the proposed toolbox indicates that the reconstruction was performed successfully. The results were evaluated using a structural similarity metric. Visual comparison against the slices obtained with the built-in commercial reconstruction software shows the potential of the proposed reconstruction pipeline, which can be easily modified to incorporate sophisticated algorithms and accelerated using parallel computation techniques.

## How to install?

 1. Download the toolbox or clone the directory [here](https://github.com/LAVI-USP/DBT-Reconstruction.git).
 
 3. The toolbox was tested on **Windows** 10 x64, **Linux** (Ubuntu) x64, and **MacoS** High Sierra machines, with Matlab 2015a and 2017b. 

Please report issues [here](https://github.com/LAVI-USP/DBT-Reconstruction/issues).

> **Note:** In order to use CUDA, make sure you have a GPU that works with CUDA and everything is updated. This [benchmark](https://www.mathworks.com/help/distcomp/examples/benchmarking-a-b-on-the-gpu.html) is very helpful to test your GPU.

## Contribute? 

We are pleased with any contributions. Fell free to make any [pull requests](https://github.com/LAVI-USP/DBT-Reconstruction/pulls) or send us an e-mail.


## Some results:

**DBT projection simulation:**

![Imgur](https://i.imgur.com/tcgqs8I.gif?1)


**Projections:**

![Imgur](https://i.imgur.com/VXk6MFi.gif?2)


**Reconstructed volume:**

![Imgur](https://i.imgur.com/R9tB35f.gif?2)


## Contact:

If you have any question or suggestion, please send us an e-mail:

- Rodrigo - rodrigo.vimieiro@gmail.com or rodrigo.vimieiro@usp.br
- Marcelo - mvieira@sc.usp.br

## License:

The toolbox is licensed under the **GNU General Public License v3.0**. Please check the [licence file](https://github.com/LAVI-USP/DBT-Reconstruction/blob/master/LICENSE).

## Reference:

If you use the toolbox, we will be very grateful if you refer to this [paper](http://iris.sel.eesc.usp.br/lavi/):

> Vimieiro, R. B., Borges, L. R., Vieira, M. A. C. (2019) Open-source reconstruction toolbox for digital
breast tomosynthesis. In: Costa-Felix R., Machado J. C., Alvarenga A. V. (eds) XXVI Brazilian Congress on Biomedical Engineering 2018. IFMBE Proceedings, vol 70/1. Springer, Singapore.

---

Laboratory of Computer Vision ([Lavi](http://iris.sel.eesc.usp.br/lavi/))  
Department of Electrical and Computer Engineering  
São Carlos School of Engineering, University of São Paulo  
São Carlos - Brazil  
