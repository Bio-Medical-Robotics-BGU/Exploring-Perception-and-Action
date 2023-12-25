''' This code uses the split into train and test data for each fold, and creates the datasets for each fold '''

# %% Imports
import os
import numpy as np
import scipy.io

# Paths
# Replace the following directory with the location of the saved features
MatPath = r"C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final"

# Replace the following directory with the location of the saved fold indices, and to which the datasets will be saved
ProjectPath = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final"
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

k = 10


# %% Loading splits
os.chdir(DatasetPath)


#load splits for the non-negative effect participants
TrainIndsSplit = np.load("Dictionary_TrainIndSplit.npy", allow_pickle = 'TRUE').item()
TestIndsSplit = np.load("Dictionary_TestIndSplit.npy", allow_pickle = 'TRUE').item()

    
 
# Taking Current Split 
for Split in np.arange(1, k + 1):
    print(Split)
    
    TrainingParticipants = TrainIndsSplit[f'split_{Split}']
    TestParticipants = TestIndsSplit[f'split_{Split}']

    # Getting Data
    os.chdir(MatPath) 
    
    #train
    for i in TrainingParticipants:
        print(i)
        #features
        os.chdir(MatPath)
        sigs = scipy.io.loadmat(f'Features_SN{i}.mat')
        sigs = sigs['AllFeatures']
        
        features_ref = sigs[:, 0, :]
        features_comp = sigs[:, 1, :]
        
      
        if i == TrainingParticipants[0]:
            TrainFeaturesComp = features_comp
            TrainFeaturesRef = features_ref
        
        else:
            TrainFeaturesComp = np.concatenate((TrainFeaturesComp, features_comp), axis = 0)
            TrainFeaturesRef = np.concatenate((TrainFeaturesRef, features_ref), axis = 0)
            
    #test
    for i in TestParticipants:
        print(i)
        #features
        os.chdir(MatPath)
        sigs = scipy.io.loadmat(f'Features_SN{i}.mat')
        sigs = sigs['AllFeatures'][:192, :, :]
        
        features_ref = sigs[:, 0, :]
        features_comp = sigs[:, 1, :]
      
    
      
        if i == TestParticipants[0]:
            TestFeaturesComp = features_comp
            TestFeaturesRef = features_ref
        
        else:
            TestFeaturesComp = np.concatenate((TestFeaturesComp, features_comp), axis = 0)
            TestFeaturesRef = np.concatenate((TestFeaturesRef, features_ref), axis = 0)
          
    os.chdir(DatasetPath)
    
    # Saving
    #train
    np.save(f'TrainFeaturesComp_Split{Split}', TrainFeaturesComp)
    np.save(f'TrainFeaturesRef_Split{Split}', TrainFeaturesRef)
    
    #test
    np.save(f'TestFeaturesComp_Split{Split}', TestFeaturesComp)
    np.save(f'TestFeaturesRef_Split{Split}', TestFeaturesRef)
 