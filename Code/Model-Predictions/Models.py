''' This code computes the responses of each of the four models'''

# %% Imports
import os
import numpy as np
import scipy.io

# %% Paths
# Replace the following directory with the project directory
base = "D:\OneDrive\PerceptionActionReview"
MatPath = os.path.join(base, "Preprocessed")

ProjectPath = os.path.join(base, "saved_predictions") 

# %% Negative
Participants = np.arange(1, 41)
Participants = np.setdiff1d(Participants, np.array([22, 40]))

for i in Participants:

  os.chdir(MatPath)
  sigs = scipy.io.loadmat(f'Features_SN{i}.mat')
  sigs = sigs['AllFeatures']
  
  features_ref = sigs[:, 0, :]
  features_comp = sigs[:, 1, :]
    
  #taking only the positive session trials - the first 192
  features_ref = features_ref[0:192, :]
  features_comp = features_comp[0:192, :]
  
  
  DistAnswers = np.zeros((len(features_ref), 1))
  MaxVelAnswers = np.zeros((len(features_ref), 1))
  MeanVelAnswers = np.zeros((len(features_ref), 1))
  MaxGFAnswers = np.zeros((len(features_ref), 1))
  
  for j in range(len(features_ref)):
    if features_ref[j, 0] < features_comp[j, 0]:
      DistAnswers[j] = 1
    
    if features_ref[j, 1] < features_comp[j, 1]:
      MaxVelAnswers[j] = 1
    
    if features_ref[j, 2] < features_comp[j, 2]:
      MeanVelAnswers[j] = 1
    
    if features_ref[j, 3] > features_comp[j, 3]:
      MaxGFAnswers[j] = 1
    

  mdic1 = {"Preds": DistAnswers, "label": "Preds"}
  mdic2 = {"Preds": MaxVelAnswers, "label": "Preds"}
  mdic3 = {"Preds": MeanVelAnswers, "label": "Preds"}
  mdic4 = {"Preds": MaxGFAnswers, "label": "Preds"}
  
  os.chdir(ProjectPath)  
  scipy.io.savemat(f'Preds_SN{i}_MaxDist.mat', mdic1)
  scipy.io.savemat(f'Preds_SN{i}_MaxVel.mat', mdic2)
  scipy.io.savemat(f'Preds_SN{i}_MeanVel.mat', mdic3)
  scipy.io.savemat(f'Preds_SN{i}_MaxGF.mat', mdic4)
