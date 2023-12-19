''' This code takes one of the runs for a network and examines the certainty values of the predictions for each of the
folds. This code creates Fig. 7 in the manuscript.'''
# %% Imports 
from tensorflow.keras.models import Model
from tensorflow.keras.layers import LSTM, Dense, Dropout, BatchNormalization, Input, Subtract, Bidirectional
from tensorflow.keras.regularizers import l2
import matplotlib.pyplot as plt
from tensorflow.keras.callbacks import LearningRateScheduler
import sklearn.metrics as sk
import tensorflow as tf
import tensorflow.keras.backend as K
from keras_self_attention import SeqSelfAttention
from tensorflow.keras.models import load_model

import os
import numpy as np
import scipy.io

plt.rcParams["font.family"] = "Times New Roman"

# %% Paths
# Replace the following directory with the location of the saved preprocessed data
MatPath = r"C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final"

# Replace the following directory with the location where the datasets are saved
ProjectPath = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final"
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

# Replace the following directory with the location where the models are saved
save_path = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final\saved_predictions"

# %% Setting up model names
run = int(input("Which run would you like to use (number between 1 - 5)? \n"))

ModelNames = {}

for f in range(1, 11):
  ModelNames[f'Fold{f}'] = f'Network_PositionVelocityAccelerationGripForce_Fold{f}_Run{run}.h5'


# %% Getting predictions for each fold
#all trials
ErrorCerts = {}
CorrectCerts = {}

#force only trials
ErrorForceCerts = {}
CorrectForceCerts = {}

#skin stretch trials
ErrorSSCerts = {}
CorrectSSCerts = {}

#running over the 10 folds
for i in range(1, 11):
  Split = i

  #getting signals
  os.chdir(DatasetPath)

  AllTestSignalsComp = np.load(f'TestSignalsComp_Split{Split}.npy')
  AllTestSignalsRef = np.load(f'TestSignalsRef_Split{Split}.npy')
  AllTestPlabels = np.load(f'TestPlabels_Split{Split}.npy')
  AllTestKComps = np.load(f'TestKComps_Split{Split}.npy')
  AllTestTdGains = np.load(f'TestTdGains_Split{Split}.npy')


  #loading model
  os.chdir(save_path)
  model_name = ModelNames[f'Fold{i}']
  model = load_model(model_name, custom_objects={'SeqSelfAttention': SeqSelfAttention})

  #getting the model outputs for the test signals
  outs = model.predict([AllTestSignalsComp, AllTestSignalsRef])
  preds = np.zeros_like(outs)
  for ii in range(outs.shape[0]):
    if outs[ii]>0.5:
      preds[ii] = 1

  #for all trials together
  correct = np.flatnonzero(preds == AllTestPlabels)
  errors = np.flatnonzero(preds != AllTestPlabels)
  
  #for force trials and skin stretch trials separately
  force_trials = np.flatnonzero(AllTestTdGains == 0)
  stretch_trials = np.flatnonzero(AllTestTdGains != 0)
    
  correct_force_trials = np.intersect1d(correct, force_trials)
  errors_force_trials = np.intersect1d(errors, force_trials)

  correct_stretch_trials = np.intersect1d(correct, stretch_trials)
  errors_stretch_trials = np.intersect1d(errors, stretch_trials)

  #check
  assert len(np.unique(np.concatenate((correct_stretch_trials, errors_stretch_trials), axis = 0))) == stretch_trials.shape[0]
  assert len(np.unique(np.concatenate((correct_stretch_trials, errors_stretch_trials, correct_force_trials,
                                       errors_force_trials), axis = 0))) == AllTestPlabels.shape[0]
  
  #splitting into the different comparison stiffness values
  kcomps = np.unique(AllTestKComps)
  
  for j, k in enumerate(kcomps):

    rel_k = np.flatnonzero(AllTestKComps == k)
    
    #all trials
    err_k = np.intersect1d(rel_k, errors)
    corr_k = np.intersect1d(rel_k, correct)
    
    err_certs_k = outs[err_k]
    corr_certs_k = outs[corr_k]
    
    if i == 1:
      ErrorCerts[f'k{k}'] = err_certs_k
      CorrectCerts[f'k{k}'] = corr_certs_k
    else:
      ErrorCerts[f'k{k}'] = np.concatenate((ErrorCerts[f'k{k}'], err_certs_k), axis = 0)
      CorrectCerts[f'k{k}'] = np.concatenate((CorrectCerts[f'k{k}'], corr_certs_k), axis = 0)
    
    
    #force trials
    err_force_k = np.intersect1d(rel_k, errors_force_trials)
    corr_force_k = np.intersect1d(rel_k, correct_force_trials)

    err_force_certs_k = outs[err_force_k]
    corr_force_certs_k = outs[corr_force_k]
    
    if i == 1:
      ErrorForceCerts[f'k{k}'] = err_force_certs_k
      CorrectForceCerts[f'k{k}'] = corr_force_certs_k
    else:
      ErrorForceCerts[f'k{k}'] = np.concatenate((ErrorForceCerts[f'k{k}'], err_force_certs_k), axis = 0)
      CorrectForceCerts[f'k{k}'] = np.concatenate((CorrectForceCerts[f'k{k}'], corr_force_certs_k), axis = 0)
    
    #stretch trials
    err_ss_k = np.intersect1d(rel_k, errors_stretch_trials)
    corr_ss_k = np.intersect1d(rel_k, correct_stretch_trials)
    
    err_ss_certs_k = outs[err_ss_k]
    corr_ss_certs_k = outs[corr_ss_k]
    
    if i == 1:
      ErrorSSCerts[f'k{k}'] = err_ss_certs_k
      CorrectSSCerts[f'k{k}'] = corr_ss_certs_k
    else:
      ErrorSSCerts[f'k{k}'] = np.concatenate((ErrorSSCerts[f'k{k}'], err_ss_certs_k), axis = 0)
      CorrectSSCerts[f'k{k}'] = np.concatenate((CorrectSSCerts[f'k{k}'], corr_ss_certs_k), axis = 0)
      
# %% Creating vectors of means for each case

# all trials
CorrectCertsVec = np.zeros((len(kcomps), 1))
ErrorCertsVec = np.zeros((len(kcomps), 1))

CorrectForceCertsVec = np.zeros((len(kcomps), 1))
ErrorForceCertsVec = np.zeros((len(kcomps), 1))

CorrectSSCertsVec = np.zeros((len(kcomps), 1))
ErrorSSCertsCertsVec = np.zeros((len(kcomps), 1))

for j, k in enumerate(kcomps):
  CorrectCertsVec[j] = np.mean(CorrectCerts[f'k{k}'])
  ErrorCertsVec[j] = np.mean(ErrorCerts[f'k{k}'])
  
  CorrectForceCertsVec[j] = np.mean(CorrectForceCerts[f'k{k}'])
  ErrorForceCertsVec[j] = np.mean(ErrorForceCerts[f'k{k}'])
  
  CorrectSSCertsVec[j] = np.mean(CorrectSSCerts[f'k{k}'])
  ErrorSSCertsCertsVec[j] = np.mean(ErrorSSCerts[f'k{k}'])
  
# %% choosing who to plot this time - FINISH
fin = input("Hanna - did you finish plots?")
CorrectCertsRel = CorrectSSCertsVec
ErrorCertsRel = ErrorSSCertsCertsVec

# %% plotting
plt.figure()
p1 = plt.scatter(0, CorrectCertsRel[0], c = 'green', marker = '*')
p2 = plt.scatter(0, ErrorCertsRel[0], c = 'red', marker = 's')

plt.scatter(np.arange(len(kcomps) - 1) + 1, CorrectCertsRel[1:], c = 'green', marker = '*')
plt.scatter(np.arange(len(kcomps) - 1) + 1, ErrorCertsRel[1:], c = 'red', marker = 's')

plt.legend([p1, p2], ['Correct', 'Errors'])

plt.xlim([-1.1, 12.1])
plt.ylim([0, 1])
plt.plot(np.arange(-1.1, 12.2, 0.1), 0.5*np.ones((len(np.arange(-1.1, 12.2, 0.1))),), linestyle='dashed', c = 'black', linewidth = 0.5)

plt.plot(5.5*np.ones((len(np.arange(0, 1.1, 0.1))),), np.arange(0, 1.1, 0.1), linestyle = ':', c = [0.2, 0.2, 0.2], linewidth = 0.5)

plt.xticks(np.arange(len(kcomps)), labels = kcomps, fontsize = 14)
plt.yticks(np.arange(0, 1.1, 0.1), fontsize = 14)

# plt.xlabel('Comparison Stiffness Level [N/m]',  fontsize = 14)
# plt.ylabel('Average Network Output',  fontsize = 14)

# %% plotting with different colors
# %% all trials
CorrectCertsRel = CorrectCertsVec
ErrorCertsRel = ErrorCertsVec

plt.figure()
p1 = plt.scatter(0, CorrectCertsRel[0], c = np.array([0, 0, 0])/255, marker = '*')
p2 = plt.scatter(0, ErrorCertsRel[0], c = np.array([166, 166, 166])/255, marker = 's')

plt.scatter(np.arange(len(kcomps) - 1) + 1, CorrectCertsRel[1:], c = np.array([0, 0, 0])/255, marker = '*')
plt.scatter(np.arange(len(kcomps) - 1) + 1, ErrorCertsRel[1:], c = np.array([166, 166, 166])/255, marker = 's')

plt.legend([p1, p2], ['Correct', 'Errors'])

plt.xlim([-1.1, 12.1])
plt.ylim([0, 1])
plt.plot(np.arange(-1.1, 12.2, 0.1), 0.5*np.ones((len(np.arange(-1.1, 12.2, 0.1))),), linestyle='dashed', c = 'black', linewidth = 0.5)

plt.plot(5.5*np.ones((len(np.arange(0, 1.1, 0.1))),), np.arange(0, 1.1, 0.1), linestyle = ':', c = [0.2, 0.2, 0.2], linewidth = 0.5)

plt.xticks(np.arange(len(kcomps)), labels = kcomps, fontsize = 14)
plt.yticks(np.arange(0, 1.1, 0.1), fontsize = 14)

# %% force trials
CorrectCertsRel = CorrectForceCertsVec
ErrorCertsRel = ErrorForceCertsVec

plt.figure()
# p1 = plt.scatter(0, CorrectCertsRel[0], c = np.array([175, 126, 5])/255, marker = '*')
# p2 = plt.scatter(0, ErrorCertsRel[0], c = np.array([249, 186, 33])/255, marker = 's')

# plt.scatter(np.arange(len(kcomps) - 1) + 1, CorrectCertsRel[1:], c = np.array([175, 126, 5])/255, marker = '*')
# plt.scatter(np.arange(len(kcomps) - 1) + 1, ErrorCertsRel[1:], c = np.array([249, 186, 33])/255, marker = 's')
p1 = plt.scatter(0, CorrectCertsRel[0], c = np.array([0, 0, 0])/255, marker = '*')
p2 = plt.scatter(0, ErrorCertsRel[0], c = np.array([166, 166, 166])/255, marker = 's')

plt.scatter(np.arange(len(kcomps) - 1) + 1, CorrectCertsRel[1:], c = np.array([0, 0, 0])/255, marker = '*')
plt.scatter(np.arange(len(kcomps) - 1) + 1, ErrorCertsRel[1:], c = np.array([166, 166, 166])/255, marker = 's')

plt.legend([p1, p2], ['Correct', 'Errors'])

plt.xlim([-1.1, 12.1])
plt.ylim([0, 1])
plt.plot(np.arange(-1.1, 12.2, 0.1), 0.5*np.ones((len(np.arange(-1.1, 12.2, 0.1))),), linestyle='dashed', c = 'black', linewidth = 0.5)

plt.plot(5.5*np.ones((len(np.arange(0, 1.1, 0.1))),), np.arange(0, 1.1, 0.1), linestyle = ':', c = [0.2, 0.2, 0.2], linewidth = 0.5)

plt.xticks(np.arange(len(kcomps)), labels = kcomps, fontsize = 14)
plt.yticks(np.arange(0, 1.1, 0.1), fontsize = 14)

# %% ss trials
CorrectCertsRel = CorrectSSCertsVec
ErrorCertsRel = ErrorSSCertsCertsVec

plt.figure()
# p1 = plt.scatter(0, CorrectCertsRel[0], c = np.array([0, 67, 112])/255, marker = '*')
# p2 = plt.scatter(0, ErrorCertsRel[0], c = np.array([0, 121, 204])/255, marker = 's')

# plt.scatter(np.arange(len(kcomps) - 1) + 1, CorrectCertsRel[1:], c = np.array([0, 67, 112])/255, marker = '*')
# plt.scatter(np.arange(len(kcomps) - 1) + 1, ErrorCertsRel[1:], c = np.array([0, 121, 204])/255, marker = 's')
p1 = plt.scatter(0, CorrectCertsRel[0], c = np.array([0, 0, 0])/255, marker = '*')
p2 = plt.scatter(0, ErrorCertsRel[0], c = np.array([166, 166, 166])/255, marker = 's')

plt.scatter(np.arange(len(kcomps) - 1) + 1, CorrectCertsRel[1:], c = np.array([0, 0, 0])/255, marker = '*')
plt.scatter(np.arange(len(kcomps) - 1) + 1, ErrorCertsRel[1:], c = np.array([166, 166, 166])/255, marker = 's')

plt.legend([p1, p2], ['Correct', 'Errors'])

plt.xlim([-1.1, 12.1])
plt.ylim([0, 1])
plt.plot(np.arange(-1.1, 12.2, 0.1), 0.5*np.ones((len(np.arange(-1.1, 12.2, 0.1))),), linestyle='dashed', c = 'black', linewidth = 0.5)

plt.plot(5.5*np.ones((len(np.arange(0, 1.1, 0.1))),), np.arange(0, 1.1, 0.1), linestyle = ':', c = [0.2, 0.2, 0.2], linewidth = 0.5)

plt.xticks(np.arange(len(kcomps)), labels = kcomps, fontsize = 14)
plt.yticks(np.arange(0, 1.1, 0.1), fontsize = 14)