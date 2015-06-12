These images come from the publicly available retinal image database DRIVE (http://www.isi.uu.nl/Research/Databases/DRIVE/) described in the paper:

J.J. Staal, M.D. Abramoff, M. Niemeijer, M.A. Viergever, B. van Ginneken, "Ridge based vessel segmentation in color images of the retina", IEEE Transactions on Medical Imaging, 2004, vol. 23, pp. 501-509.

The 20 test images of the base have been preprocessed using the following procedures:
1) Extract the green channel;
2) Perform a black top-hat by a circular structuring element of 5 pixels diameter;
3) Set all pixels outside the mask of the eye fundus to 0.

