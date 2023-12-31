''' This is the neural network code using only position and grip force, and performs the ablation analysis'''
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

# Replace the following directory with the location of the saved preprocessed data
MatPath = r"D:\OneDrive\PerceptionAction\Preprocessed"

# Replace the following directory with the location where the datasets are saved
ProjectPath = r"D:\OneDrive\PerceptionAction\Data"
DatasetPath = os.path.join(ProjectPath, 'DatasetsFolds')

# Replace the following directory with the location where the model predictions should be saved
save_path = r"D:\OneDrive\PerceptionAction\saved_predictions"

Splits = np.arange(1, 11) #for the k-fold validation
Runs = np.arange(1, 6) #for the 5 repititions


# %% runs neural network
#include only position and grip force
sigs = np.array([0, 3])

abl = int(input("Which layers should be included? \n 1. Block #1 \n 2. Block #1 + last LSTM \n 3. Block #1 + attention + last LSTM \n 4. Full model \n"))

save_model = int(input("Would you like to save the trained networks? \n 1. Yes \n 2. No \n"))

for run in Runs:
  for Split in Splits:
    os.chdir(DatasetPath)
    
    #train
    AllTrainSignalsComp = np.load(f'TrainSignalsComp_Split{Split}.npy')
    AllTrainSignalsRef = np.load(f'TrainSignalsRef_Split{Split}.npy')
    AllTrainPlabels = np.load(f'TrainPlabels_Split{Split}.npy')

    #Test
    AllTestSignalsComp = np.load(f'TestSignalsComp_Split{Split}.npy')
    AllTestSignalsRef = np.load(f'TestSignalsRef_Split{Split}.npy')
    AllTestPlabels = np.load(f'TestPlabels_Split{Split}.npy')
    
    # taking only position and grip force
    AllTrainSignalsComp = AllTrainSignalsComp[:, :, sigs]
    AllTrainSignalsRef = AllTrainSignalsRef[:, :, sigs]
    
    AllTestSignalsComp = AllTestSignalsComp[:, :, sigs]
    AllTestSignalsRef = AllTestSignalsRef[:, :, sigs]
    
    
    # Building neural network
    lstm1 = 128
    lstm2 = 64 
    lstm3 = 64

    do = 0.4 
    reg = 0.001

    K.clear_session()
    
    Input1 = Input(shape = (AllTrainSignalsComp.shape[1], AllTrainSignalsComp.shape[2]))
    Input2 = Input(shape = (AllTrainSignalsRef.shape[1], AllTrainSignalsRef.shape[2]))
    
    #only block 1
    if abl == 1: #only block 1
      x1 = Bidirectional(LSTM(lstm1, kernel_regularizer=l2(reg), return_sequences = False))(Input1)
      x2 = Dropout(do)(x1)
      x3 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(x2)
      
      y1 = Bidirectional(LSTM(lstm1, kernel_regularizer=l2(reg), return_sequences = False))(Input2)
      y2 = Dropout(do)(y1)
      y3 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(y2)
      
      subtracted = Subtract(name = 'final_rep')([y3, x3])
      
    else: #options 2 - 4
      #included in options 2 - 4:
      x1 = Bidirectional(LSTM(lstm1, kernel_regularizer=l2(reg), return_sequences = True))(Input1)
      x2 = Dropout(do)(x1)
      x3 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(x2)
      
      y1 = Bidirectional(LSTM(lstm1, kernel_regularizer=l2(reg), return_sequences = True))(Input2)
      y2 = Dropout(do)(y1)
      y3 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(y2)
    
      if abl == 2: #block 1 + last LSTM
        x8 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='CompRep')(x3)
        y8 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='RefRep')(y3)
        
      elif abl == 3: #block 1 + attention + last LSTM
        x4 = SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL,
                        attention_activation='softmax',
                        name='Attention1')(x3)
        x8 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='CompRep')(x4)
        
        y4 = SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL,
                        attention_activation='softmax',
                        name='Attention2')(y3)
        y8 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='RefRep')(y4)
        
      else: #full model
        
        x4 = SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL,
                        attention_activation='softmax',
                        name='Attention1')(x3)
        x5 = Bidirectional(LSTM(lstm2, return_sequences = True, kernel_regularizer=l2(reg)))(x4)
        x6 = Dropout(do)(x5)
        x7 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(x6)
        x8 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='CompRep')(x7)
    
        
        y4 = SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL,
                        attention_activation='softmax',
                        name='Attention2')(y3)
        y5 = Bidirectional(LSTM(lstm2, return_sequences = True, kernel_regularizer=l2(reg)))(y4)
        y6 = Dropout(do)(y5)
        y7 = BatchNormalization(axis = -1, momentum = 0.99, epsilon = 0.001, center = True, scale = True)(y6)
        y8 = Bidirectional(LSTM(lstm3, return_sequences = False, kernel_regularizer=l2(reg)), name='RefRep')(y7)

      subtracted = Subtract(name = 'final_rep')([y8, x8])
  
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
      
    all_abls = ['Block1', 'Block1_LastLSTM', 'Block1_Att_LastLSTM', 'FullModel']


    if save_model == 1:
      # Replace the following directory with the location where the model should be saved
      os.chdir(save_path)  
      model_name = f'Network_{all_abls[abl - 1]}_Split{Split}_Run{run}.h5'
      model.save(model_name)
      
    # getting model predictions for testing participants of this fold
    os.chdir(DatasetPath)
    TestIndsSplit = np.load("Dictionary_TestIndSplit.npy", allow_pickle = 'TRUE').item()
    
    TestParticipants = TestIndsSplit[f'split_{Split}']

    for i in TestParticipants:
      os.chdir(MatPath)

      #comparison
      comps = scipy.io.loadmat(f'CompSignals_SN{i}.mat')
      comps = comps['AllCompSigs']
      
      #standard
      refs = scipy.io.loadmat(f'StandardSignals_SN{i}.mat')
      refs = refs['AllRefSigs']
     
      comps = comps[:192, :, sigs]
      refs = refs[:192, :, sigs]
      
      outs = model.predict([comps, refs])
      preds = np.zeros_like(outs)
      for j in range(outs.shape[0]):
        if outs[j]>0.5:
          preds[j] = 1
      
      mdic1 = {"Preds": preds, "label": "Preds"}
      
      #saving predictions      
      os.chdir(save_path)  
      
      scipy.io.savemat(f'Preds_SN{i}_{all_abls[abl - 1]}_Run{run}.mat', mdic1)





