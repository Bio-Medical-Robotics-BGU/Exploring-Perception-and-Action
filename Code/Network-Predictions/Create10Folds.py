''' This code creates the participant folds for k-fold validation - that is, defines which participants will be the train and 
test participants for each of the 10 folds'''

# %% Imports
from sklearn.model_selection import KFold
import os
import numpy as np
from scipy.io import savemat

# %% Paths
# Replace the following directory with the location of the saved preprocessed data
MatPath = r"C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final"

# Replace the following directory with the location to which to save the fold indices
ProjectPath = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final"
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

k = 10
# %% Without the four negative effect participants
Participants = np.arange(1, 41)
Participants = np.setdiff1d(Participants, np.array([2, 4, 5, 21, 22, 40])) 
#participants 22 and 40 were the outlier participants that were removed from all analyses
#participants 2, 4, 5, 21 were the negative effect participants

kf = KFold(n_splits = k, shuffle = True)

TrainIndsSplit = {}
TestIndsSplit = {}

i = 1
for train_index, test_index in kf.split(Participants):
 
  TrainIndsSplit[f'split_{i}'] = Participants[train_index]
  TestIndsSplit[f'split_{i}'] = Participants[test_index]

  i += 1
  
os.chdir(DatasetPath)
np.save("Dictionary_TrainIndSplit.npy", TrainIndsSplit)
np.save("Dictionary_TestIndSplit.npy", TestIndsSplit)

#to save for opening in MATLAB
savemat('TestIndsSplit.mat', TestIndsSplit, oned_as='row')

# %% Randomly adding in the four negative effect participants
AddParticipants = np.array([2, 4, 5, 21])

#adding them into the folds while ensuring same size folds as much as possible
Fours = np.array([1, 2, 3, 4]) #the folds with three participants in the test data (30 in the train data)
Threes = np.array([5, 6, 7, 8, 9, 10]) #the folds with four participants in the test data (31 in the train data))

for i in Fours:
  TrainIndsSplit[f'split_{i}'] = np.concatenate((TrainIndsSplit[f'split_{i}'], AddParticipants))
  
#adding each of the neg effect participants to one of the Three folds randomly:
AddTo = np.random.choice(Threes, len(AddParticipants), replace = False)
AddTo = np.sort(AddTo)
count = 0

for i in Threes:
  if i == AddTo[count]:
    TestIndsSplit[f'split_{i}'] = np.concatenate((TestIndsSplit[f'split_{i}'], np.reshape(np.array(AddParticipants[count]), (1,))))
    TrainIndsSplit[f'split_{i}'] = np.concatenate((TrainIndsSplit[f'split_{i}'], np.setdiff1d(AddParticipants, AddParticipants[count])))
    count += 1
  else:
    TrainIndsSplit[f'split_{i}'] = np.concatenate((TrainIndsSplit[f'split_{i}'], AddParticipants))
    

os.chdir(DatasetPath)
np.save("Dictionary_TrainIndSplit_AllParticipants.npy", TrainIndsSplit)
np.save("Dictionary_TestIndSplit_AllParticipants.npy", TestIndsSplit)

#to save for opening in MATLAB
savemat('TestIndsSplit_AllParticipants.mat', TestIndsSplit, oned_as='row')