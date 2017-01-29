# A Cloud Intrusion Detection System Using Novel PRFCM Clustering and KNN Based Dempster-Shafer Rule

DOI: [10.4018/IJCAC.2016100102](http://www.igi-global.com/article/a-cloud-intrusion-detection-system-using-novel-prfcm-clustering-and-knn-based-dempster-shafer-rule/173770)

An Anomaly based method for implementing Intrusion Detection Systems is described which uses a Novel PRFCM Clustering and KNN based Dempster-Shafer Rule. The PRFCM Clustering is a modification to the existing FCM Clustering Algortihm. Pre-Processed NSL-KDD dataset was used to get the results. The results are evaluated and compared with original FCM algortihm in the paper.

## Introduction

The FCM algorithm is very sensitive to noise and is easily struck at local optima. In order to overcome
this limitation, spatial context of connection is taken into account considering its neighboring
connections. A Penalty Reward based FCM algorithm is implemented here which can handle small
as well as large amount of noise by adjusting a penalty and reward coefficient. The algorithm takes
into account both the feature information and spatial information. The objective function is
modified to incorporate the penalty and reward term by which it can overcome the local optima. The new objective function of the PRFCM algorithm is defined as follows:

![alt tag](https://github.com/shakti365/IDS-PRFCM-Clustering/blob/master/resources/PRFCM_objfun.png)

### Dataset

The experiments are performed on NSL-KDD Train and Test Dataset. These dataset were pre-processed and normalized before use. It can be obtained from the following source. 

[NSL-KDD Dataset](http://www.unb.ca/research/iscx/dataset/iscx-NSL-KDD-dataset.html)

### Results

The following results show the performance of PRFCM clustering over FCM clustering algorithm

![alt tag](https://github.com/shakti365/IDS-PRFCM-Clustering/blob/master/resources/f2d535ba8c6818c945eede519011e995-original.jpg) 

## Authors

* **Partha Ghosh**
* **Shivam Shakti**
* **Santanu Phadikar**

## Copyright

This paper was published in [International Journal of Cloud Applications and Computing](http://www.igi-global.com/journal/international-journal-cloud-applications-computing/41974)

Volume 6 • Issue 4 • October-December 2016

