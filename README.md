# EEWCAN

# Collaborative Attention Network for Rapid Epicentral Distance and Magnitude Estimation
Jingbao Zhu, Shanyou Li, Qiang Ma, Kunpeng Yao, Pengjie Huang, and Jindong Song



# Data
Due to the large memory size of waveform data and limitations in data transmission, 
we uploaded the waveform dataset to IEEE Dataport. The datasets used for training and testing EEWCAN model 
and baseline models(CNN-Dis、CRNN-Dis、CNN-Mag、LSTM-Mag、MEANet-Mag、MagNet、EqViT), 
as well as Chinese strong-ground motion test dataset, can be downloaded from IEEE Dataport. Here is the link to IEEE Dataport:
https://ieee-dataport.org/documents/dataset-training-and-testing-eewcan.
If the Dataset help your research, please kindly cite:
Jingbao Zhu, Jindong Song, Shanyou Li (2025). Dataset for training and testing EEWCAN. IEEE Dataport. https://dx.doi.org/10.21227/wr8h-ts82.

# Data acquisition and preprocessing
These codes are used for acquisition and preprocessing of waveform data recorded by K-NET stations, in order to establish the dataset for this study.
The waveform data of K-NET network can be downloaded from this website：https://www.kyoshin.bosai.go.jp/kyoshin/quake/index_en.html

# Baseline models and EEWVAN
We provide the training and inference scripts for the baseline models and EEWCAN, the training codes for these models, 
and the pre training weights for each trained baseline models and EEWCAN

# requirements.txt
provided environment or dependency specifications

# evaluation.ipynb
We provide codes for model performance evaluation
