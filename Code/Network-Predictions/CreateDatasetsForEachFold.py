''' This code uses the split into train and test data for each fold, and creates the datasets for each fold '''

# %% Imports
import os
import numpy as np
import scipy.io

# Paths
# Replace the following directory with the location of the saved preprocessed data
MatPath = r"C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final"

# Replace the following directory with the location of the saved fold indices, and to which the datasets will be saved
ProjectPath = r"C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final"
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

k = 10

# %% Loading splits
os.chdir(DatasetPath)

flag = 1
Dataset = int(input("Which dataset would you like to create? \n 1. Dataset 1 (non-negative participants) \n 2. Dataset 2 (all participants) \n 3. Dataset 3 (both stretch conditions)"))

while (flag):
    if Dataset == 1:
        #load splits for the non-negative effect participants
        TrainIndsSplit = np.load("Dictionary_TrainIndSplit.npy", allow_pickle = 'TRUE').item()
        TestIndsSplit = np.load("Dictionary_TestIndSplit.npy", allow_pickle = 'TRUE').item()
        flag = 0
    
    elif Dataset == 2 or Dataset == 3:
        #load splits for all participants (this split is also used for Dataset 3, which includes both the positive and 
        #negative stretch trials for all the participants)
        TrainIndsSplit = np.load("Dictionary_TrainIndSplit_AllParticipants.npy", allow_pickle = 'TRUE').item()
        TestIndsSplit = np.load("Dictionary_TestIndSplit_AllParticipants.npy", allow_pickle = 'TRUE').item()
        flag = 0
    
    else:
        print("Please input 1, 2 or 3 \n")
        Dataset = int(input("Which dataset would you like to create? \n 1. Dataset 1 (non-negative participants) \n 2. Dataset 2 (all participants) \n 3. Dataset 3 (both stretch conditions)"))
        flag = 1
        
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
        #signals
        if Dataset == 1 or Dataset == 2:
            sigs_comp = scipy.io.loadmat(f'CompSignals_SN{i}.mat') 
        else:
            #for both positive and negative stretch:
            sigs_comp = scipy.io.loadmat(f'CompSignals_BothStretch_SN{i}.mat') 
        sigs_comp = sigs_comp['AllCompSigs']
      
        if Dataset == 1 or Dataset == 2:
            sigs_ref = scipy.io.loadmat(f'StandardSignals_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            sigs_ref = scipy.io.loadmat(f'StandardSignals_BothStretch_SN{i}.mat') 
        sigs_ref = sigs_ref['AllRefSigs']
      
        #participant answers
        if Dataset == 1 or Dataset == 2:
            panswers = scipy.io.loadmat(f'Labels_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            panswers = scipy.io.loadmat(f'Labels_BothStretch_SN{i}.mat') 
        panswers = panswers['AllPlabels']
      
        #Kcomps
        if Dataset == 1 or Dataset == 2:
            kcomps = scipy.io.loadmat(f'Kcomps_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            kcomps = scipy.io.loadmat(f'Kcomps_BothStretch_SN{i}.mat') 
        kcomps = kcomps['AllKcomps']
      
        #Tactor discplacement gains
        if Dataset == 1 or Dataset == 2:
            tdgains = scipy.io.loadmat(f'TdGains_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            tdgains = scipy.io.loadmat(f'TdGains_BothStretch_SN{i}.mat') 
        tdgains = tdgains['AllTdGains']
      
        assert (sigs_comp.shape[0] == sigs_ref.shape[0] == len(panswers) == len(kcomps) == len(tdgains))
                      
      
        if i == TrainingParticipants[0]:
            TrainSigsCompNeg = sigs_comp
            TrainSigsRefNeg = sigs_ref
            TrainPlabelsNeg = panswers
            TrainKCompsNeg = kcomps
            TrainTdGainsNeg = tdgains
        
        else:
            TrainSigsCompNeg = np.concatenate((TrainSigsCompNeg, sigs_comp), axis = 0)
            TrainSigsRefNeg = np.concatenate((TrainSigsRefNeg, sigs_ref), axis = 0)
            TrainPlabelsNeg = np.concatenate((TrainPlabelsNeg, panswers), axis = 0)
            TrainKCompsNeg = np.concatenate((TrainKCompsNeg, kcomps), axis = 0)
            TrainTdGainsNeg = np.concatenate((TrainTdGainsNeg, tdgains), axis = 0)
    
    #test
    for i in TestParticipants:
        print(i)
        #signals
        if Dataset == 1 or Dataset == 2:
            sigs_comp = scipy.io.loadmat(f'CompSignals_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            sigs_comp = scipy.io.loadmat(f'CompSignals_BothStretch_SN{i}.mat') 
        sigs_comp = sigs_comp['AllCompSigs']
      
        if Dataset == 1 or Dataset == 2:
            sigs_ref = scipy.io.loadmat(f'StandardSignals_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            sigs_ref = scipy.io.loadmat(f'StandardSignals_BothStretch_SN{i}.mat') 
        sigs_ref = sigs_ref['AllRefSigs']
      
        #participant answers
        if Dataset == 1 or Dataset == 2:
            panswers = scipy.io.loadmat(f'Labels_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            panswers = scipy.io.loadmat(f'Labels_BothStretch_SN{i}.mat') 
        panswers = panswers['AllPlabels']
      
        #Kcomps
        if Dataset == 1 or Dataset == 2:
            kcomps = scipy.io.loadmat(f'Kcomps_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            kcomps = scipy.io.loadmat(f'Kcomps_BothStretch_SN{i}.mat') 
        kcomps = kcomps['AllKcomps']
      
        #Tactor discplacement gains
        if Dataset == 1 or Dataset == 2:
            tdgains = scipy.io.loadmat(f'TdGains_SN{i}.mat')
        else:
            #for both positive and negative stretch:
            tdgains = scipy.io.loadmat(f'TdGains_BothStretch_SN{i}.mat') 
        tdgains = tdgains['AllTdGains']
      
        assert (sigs_comp.shape[0] == sigs_ref.shape[0] == len(panswers) == len(kcomps) == len(tdgains))
      
    
      
        if i == TestParticipants[0]:
            TestSigsCompNeg = sigs_comp
            TestSigsRefNeg = sigs_ref
            TestPlabelsNeg = panswers
            TestKCompsNeg = kcomps
            TestTdGainsNeg = tdgains
        
        else:
            TestSigsCompNeg = np.concatenate((TestSigsCompNeg, sigs_comp), axis = 0)
            TestSigsRefNeg = np.concatenate((TestSigsRefNeg, sigs_ref), axis = 0)
            TestPlabelsNeg = np.concatenate((TestPlabelsNeg, panswers), axis = 0)
            TestKCompsNeg = np.concatenate((TestKCompsNeg, kcomps), axis = 0)
            TestTdGainsNeg = np.concatenate((TestTdGainsNeg, tdgains), axis = 0)
    
    
    os.chdir(DatasetPath)
    
    if Dataset == 1:
        # Saving - for non-negative participants 
        #train
        np.save(f'TrainSignalsComp_Split{Split}', TrainSigsCompNeg)
        np.save(f'TrainSignalsRef_Split{Split}', TrainSigsRefNeg)
        
        np.save(f'TrainPlabels_Split{Split}', TrainPlabelsNeg)
        np.save(f'TrainKComps_Split{Split}', TrainKCompsNeg)
        np.save(f'TrainTdGains_Split{Split}', TrainTdGainsNeg)
        
        #test
        np.save(f'TestSignalsComp_Split{Split}', TestSigsCompNeg)
        np.save(f'TestSignalsRef_Split{Split}', TestSigsRefNeg)
        
        np.save(f'TestPlabels_Split{Split}', TestPlabelsNeg)
        np.save(f'TestKComps_Split{Split}', TestKCompsNeg)
        np.save(f'TestTdGains_Split{Split}', TestTdGainsNeg)
    
    if Dataset == 2:
        # Saving - for all participants
        #train
        np.save(f'TrainSignalsComp_AllParticipants_Split{Split}', TrainSigsCompNeg)
        np.save(f'TrainSignalsRef_AllParticipants_Split{Split}', TrainSigsRefNeg)
        
        np.save(f'TrainPlabels_AllParticipants_Split{Split}', TrainPlabelsNeg)
        np.save(f'TrainKComps_AllParticipants_Split{Split}', TrainKCompsNeg)
        np.save(f'TrainTdGains_AllParticipants_Split{Split}', TrainTdGainsNeg)
        
        #test
        np.save(f'TestSignalsComp_AllParticipants_Split{Split}', TestSigsCompNeg)
        np.save(f'TestSignalsRef_AllParticipants_Split{Split}', TestSigsRefNeg)
        
        np.save(f'TestPlabels_AllParticipants_Split{Split}', TestPlabelsNeg)
        np.save(f'TestKComps_AllParticipants_Split{Split}', TestKCompsNeg)
        np.save(f'TestTdGains_AllParticipants_Split{Split}', TestTdGainsNeg)
        
    
    if Dataset == 3:
        # Saving - for all participants with both positive and negative stretch
        #train
        np.save(f'TrainSignalsComp_BothStretch_Split{Split}', TrainSigsCompNeg)
        np.save(f'TrainSignalsRef_BothStretch_Split{Split}', TrainSigsRefNeg)
        
        np.save(f'TrainPlabels_BothStretch_Split{Split}', TrainPlabelsNeg)
        np.save(f'TrainKComps_BothStretch_Split{Split}', TrainKCompsNeg)
        np.save(f'TrainTdGains_BothStretch_Split{Split}', TrainTdGainsNeg)
        
        #test
        np.save(f'TestSignalsComp_BothStretch_Split{Split}', TestSigsCompNeg)
        np.save(f'TestSignalsRef_BothStretch_Split{Split}', TestSigsRefNeg)
        
        np.save(f'TestPlabels_BothStretch_Split{Split}', TestPlabelsNeg)
        np.save(f'TestKComps_BothStretch_Split{Split}', TestKCompsNeg)
        np.save(f'TestTdGains_BothStretch_Split{Split}', TestTdGainsNeg)
