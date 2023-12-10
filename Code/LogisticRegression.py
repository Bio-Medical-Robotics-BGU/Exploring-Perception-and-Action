''' This code computes the linear regression of the four model features'''
# %% Imports 
import sklearn.metrics as sk
from sklearn.linear_model import LogisticRegression

import os
import numpy as np
import scipy.io

# %% Paths

# Replace the following directory with the location of the saved features
MatPath = r"C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final"

# Replace the following directory with the location where the datasets are saved
ProjectPath = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final"
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

# Replace the following directory with the location where the model predictions should be saved
save_path = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final\saved_predictions"

Splits = np.arange(1, 11)

# %%
for Split in Splits:
  os.chdir(DatasetPath)
  
  #train
  TrainFeaturesComp = np.load(f'TrainFeaturesComp_Split{Split}.npy')
  TrainFeaturesRef = np.load(f'TrainFeaturesRef_Split{Split}.npy')
  TrainPlabels = np.load(f'TrainPlabels_Split{Split}.npy')

  #Test
  TestFeaturesComp = np.load(f'TestFeaturesComp_Split{Split}.npy')
  TestFeaturesRef = np.load(f'TestFeaturesRef_Split{Split}.npy')
  TestPlabels = np.load(f'TestPlabels_Split{Split}.npy')
  
  SubtractedTrain = TrainFeaturesRef - TrainFeaturesComp
  TrainPlabels = np.reshape(TrainPlabels, (TrainPlabels.shape[0],))
  
  SubtractedTest = TestFeaturesRef - TestFeaturesComp
  TestPlabels = np.reshape(TestPlabels, (TestPlabels.shape[0],))
 

  logisticRegr = LogisticRegression()
  logisticRegr.fit(SubtractedTrain, TrainPlabels)
  
  # getting model predictions for testing participants of this fold
  os.chdir(DatasetPath)
  TestIndsSplit = np.load("Dictionary_TestIndSplit.npy", allow_pickle = 'TRUE').item()

  TestParticipants = TestIndsSplit[f'split_{Split}']

  #validation

  for i in TestParticipants:
    
    #features
    os.chdir(MatPath)
    sigs = scipy.io.loadmat(f'Features_SN{i}.mat')
    sigs = sigs['AllFeatures']
    
    features_ref = sigs[:, 0, :]
    features_comp = sigs[:, 1, :]
        
    
    subbed = features_ref - features_comp 
    preds = logisticRegr.predict(subbed)
    

    mdic1 = {"Preds": preds, "label": "Preds"}

    os.chdir(save_path)  
    scipy.io.savemat(f'Preds_SN{i}_LogisticRegression.mat', mdic1)