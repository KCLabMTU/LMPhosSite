# LMPhosSite
A deep learning-based approach for general protein phosphorylation site prediction using embeddings from local window sequence and pre-trained Protein Language Model 


## Evaluate Model
Programs were executed using anaconda version: 2020.07, recommended to install the same

The programs were developed in the following environment. python : 3.8.3.final.0, python-bits : 64, OS : Linux, OS-release : 5.8.0-38-generic, machine : x86_64, processor : x86_64, pandas : 1.0.5, numpy : 1.18.5, pip : 20.1.1, scipy : 1.4.1, scikit-learn : 0.23.1., keras : 2.4.3, tensorflow : 2.3.1.

1. Please place the ST_Independent_Testing.ipynb, DeepPSP_81_window_Test_Positive_18553 (1).fasta, DeepPSP_81_window_Test_Negative_101944 (2).fasta, Index_DeepPSP_33_window_Test_Negative.txt, Index_DeepPSP_33_window_Test_Positive.txt, and Again_job_was_performed_stacked_generilization_model_5112305.h5 in the same directory where you will execute the python program. Execute the ST_Independent_Testing.ipynb python program to obtain the reported result for ST phosphorylation sites.

2. Please place the Y_Independent_Testing.ipynb, test_positive_Y_81_window_3248.fasta, test_negative_Y_81_window_14503.fasta, 14503_Index_DeepPSP_33_window_Y___Test_Negative.txt, 3248_Index_DeepPSP_33_window_Y___Test_Positive.txt, and Y_______Again_job_was_performed_stacked_generilization_model_2899160.h5 in the same directory where you will execute the python program. Execute the Y_Independent_Testing.ipynb python program to obtain the reported result for Y phosphorylation sites.


*** For your convenience we have uploaded the ProtT5 feature file extraction program for each proteins (Prot_T5_file_generation_of_Proteins.py) which will be driven by shell script (Prot_T5_file_generation_of_Proteins.sh). We have uploaded ProtT5 two sample files of two proteins in LMPhosSite/data/ProtT5 Example folder for your convenience. 

*** As well as we have uploaded corresponding 1024 feature vector extraction program for phosphorylation residue/token ('S,T,Y') (DeepPSP Positive ST feature extraction Training.ipynb) from the ProtT5 file. ***

*** We have used BeoShock High-Performance Computing resources (200 GB RAM, 18 CPUs, and 1 Nvidia Tesla V100 16GB GPU)  for ProtT5 file generation of each proteins so we suggest to use similar configuration/system to extract ProtT5 feature file. 

*** All preliminary data requred are uploaded at https://drive.google.com/drive/u/1/folders/1KD69lBFJ0qtDojza1pwyJco6rniizsD6

If you need any futher help please contact Dr. Dukka B. KC at dbkc@mtu.edu.




