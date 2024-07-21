''' This is the neural network code using only position and grip force, and includes the negative effect participants with 
or without the negative stretch session'''
# %% Imports 
from tensorflow.keras.models import Model
from tensorflow.keras.layers import LSTM, Dense, Dropout, BatchNormalization, Input, Subtract, Bidirectional
from tensorflow.keras.regularizers import l2
from tensorflow.keras.callbacks import LearningRateScheduler
import tensorflow as tf
import tensorflow.keras.backend as K
from keras_self_attention import SeqSelfAttention

import os
import numpy as np
import scipy.io

# %% Paths

# Replace the following directory with the project directory
base = "D:\OneDrive\PerceptionActionReview"
MatPath = os.path.join(base, "Preprocessed")

# the location where the datasets are saved
ProjectPath = os.path.join(base, "Data") 
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

# the location where the model predictions should be saved
save_path = os.path.join(base, "saved_predictions") 
if not os.path.exists(save_path):
    os.makedirs(save_path)
  
Splits = np.arange(1, 11) #for the k-fold validation
Runs = np.arange(1, 6) #for the 5 repititions


# %% runs neural network
#include only position and grip force
sigs = np.array([0, 3])

save_model = int(input("Would you like to save the trained networks? \n 1. Yes \n 2. No \n"))


for run in Runs:
  for Split in Splits:
    os.chdir(DatasetPath)
    
    
    #train
    AllTrainSignalsComp = np.load(f'TrainSignalsComp_BothStretch_Split{Split}.npy')
    AllTrainSignalsRef = np.load(f'TrainSignalsRef_BothStretch_Split{Split}.npy')
    AllTrainPlabels = np.load(f'TrainPlabels_BothStretch_Split{Split}.npy')

    #Test
    AllTestSignalsComp = np.load(f'TestSignalsComp_BothStretch_Split{Split}.npy')
    AllTestSignalsRef = np.load(f'TestSignalsRef_BothStretch_Split{Split}.npy')
    AllTestPlabels = np.load(f'TestPlabels_BothStretch_Split{Split}.npy')

    
    # taking the desired combination of signals
    AllTrainSignalsComp = AllTrainSignalsComp[:, :, sigs]
    AllTrainSignalsRef = AllTrainSignalsRef[:, :, sigs]
    
    AllTestSignalsComp = AllTestSignalsComp[:, :, sigs]
    AllTestSignalsRef = AllTestSignalsRef[:, :, sigs]
    
    
    # Building neural network
    lstm1 = 128
    lstm3 = 64

    do = 0.4 
    reg = 0.001

    K.clear_session()

    Input1 = Input(shape = (AllTrainSignalsComp.shape[1], AllTrainSignalsComp.shape[2]))
    x1 = Bidirectional(LSTM(lstm1, kernel_regularizer=l2(reg), return_sequences = True))(Input1)
    x2 = Dropout(do)(x1)
    x3 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(x2)


    x4 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='CompRep')(x3)

    Input2 = Input(shape = (AllTrainSignalsRef.shape[1], AllTrainSignalsRef.shape[2]))
    y1 = Bidirectional(LSTM(lstm1, kernel_regularizer=l2(reg), return_sequences = True))(Input2)
    y2 = Dropout(do)(y1)
    y3 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(y2)


    y4 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='RefRep')(y3)

    subtracted = Subtract(name = 'final_rep')([y4, x4])

    out = Dense(1, activation='sigmoid')(subtracted)
      
    model = Model(inputs = [Input1, Input2], outputs = out)
         
        

    def lr_scheduler(epoch, lr):
      return lr * 0.95


    opt = tf.keras.optimizers.Adam(learning_rate = 0.0001) 

    callbacks = [LearningRateScheduler(lr_scheduler, verbose=1)]

    model.compile(loss = 'binary_crossentropy', optimizer = opt, metrics = ['acc'])

    model.summary()

    history = model.fit([AllTrainSignalsComp, AllTrainSignalsRef], AllTrainPlabels, 
                        validation_data = ([AllTestSignalsComp, AllTestSignalsRef], AllTestPlabels), 
                        batch_size = 512, 
                        shuffle = True, epochs = 50, verbose = 1, callbacks = callbacks)
   
    #to save model:
      
    if save_model == 1:
      # Replace the following directory with the location where the model predictions should be saved
      os.chdir(save_path)  
      model_name = f'Network_BothStretch_Split{Split}_Run{run}.h5'
      model.save(model_name)
      
    # getting model predictions for testing participants of this fold
    os.chdir(ProjectPath)
    TestIndsSplit = np.load("Dictionary_TestIndSplit_AllParticipants.npy", allow_pickle = 'TRUE').item()
    
    TestParticipants = TestIndsSplit[f'split_{Split}']

    for i in TestParticipants:
      os.chdir(MatPath)
      
      #comparison
      comps = scipy.io.loadmat(f'CompSignals_BothStretch_SN{i}.mat')
      comps = comps['AllCompSigs']
        
      #standard
      refs = scipy.io.loadmat(f'StandardSignals_BothStretch_SN{i}.mat')
      refs = refs['AllRefSigs']
    
      comps = comps[:, :, sigs]
      refs = refs[:, :, sigs]
      
      outs = model.predict([comps, refs])
      preds = np.zeros_like(outs)
      for j in range(outs.shape[0]):
        if outs[j] > 0.5:
          preds[j] = 1
      
      mdic1 = {"Preds": preds, "label": "Preds"}
      
      #saving predictions      
      os.chdir(save_path)  
      
      scipy.io.savemat(f'Preds_SN{i}_BothStretch_Run{run}.mat', mdic1)




