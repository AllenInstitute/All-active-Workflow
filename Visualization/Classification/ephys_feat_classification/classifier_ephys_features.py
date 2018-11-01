#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 15:22:01 2018

@author: anin
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import json
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from ggplot import *


cwd = os.getcwd()
df = pd.DataFrame()


def get_ephys_features_df(features, protocols, feature_set, cell_metadata, df):
    
    # delete the DB_check stimulus
    
    if 'DB_check_DC' in protocols.keys():
        del protocols['DB_check_DC']
    temp_protocols = {trace_name:values['stimuli'][0]['amp'] for trace_name,values in protocols.items()}
    
    # sort the traces according to stimulus amp and ignore the non spiking trace
    sorted_traces = sorted(temp_protocols, key=temp_protocols.__getitem__)[-3:] 
    
    feature_set_revised = sorted([feat for feat in feature_set['somatic_features'] if feat not in  
                                  ['depol_block', 'check_AISInitiation' ]])
    
    ephys_feat_list = list()
    training_feat_names = list()
    None_val = False
    
   
        
    
    for i,trace in enumerate(sorted_traces):
        if i == 0:
           feature = 'mean_frequency'
           ephys_feat_list.append(features[trace]['soma'][feature][0])
           training_feat_names.append(feature+'_'+str(i))
           continue
        for feature in feature_set_revised:
            if feature in features[trace]['soma'].keys():
                ephys_feat_list.append(features[trace]['soma'][feature][0])
                training_feat_names.append(feature+'_'+str(i))
            else:
                ephys_feat_list.append(None)
                None_val = None_val or True
    
    label_feature_name = ['Cre_line']
    IN_Crelines = ['Pvalb', 'Htr3a', 'Sst']
    cell_Creline = cell_metadata['Cre_line']
    
    for IN_Creline in IN_Crelines:
        if cell_Creline.startswith(IN_Creline):
             cell_Creline = IN_Creline
             ephys_feat_list.append(cell_Creline)
             break
    
    if cell_Creline not in IN_Crelines:
        ephys_feat_list.append('Pyr')
    
   
    
    if not None_val:
        if df.shape[0] == 0:
            df=pd.DataFrame(columns = training_feat_names + label_feature_name)
            df.loc[0] = ephys_feat_list
        else:
            df.loc[df.shape[0]+1] = ephys_feat_list
    
        
    return df, training_feat_names, label_feature_name


for subdir, dirs, files in os.walk(cwd):
    if 'cell_metadata.json' in files:
        for file_ in files:
            filepath =  os.path.join(subdir, file_)
            if file_.startswith('protocols'):
                protocols = json.load(open(filepath, 'r'))
            elif file_.startswith('features'):
                features = json.load(open(filepath, 'r'))
            elif file_.startswith('feature_set'):
                feature_set = json.load(open(filepath, 'r'))
            elif file_.startswith('cell_metadata'):
                cell_metadata = json.load(open(filepath, 'r'))
    
        df, training_feat_names, label_feature_name = get_ephys_features_df(features, protocols, feature_set, cell_metadata, df)

all_cre_lines = list(set(df['Cre_line'].tolist()))
    
tsne = TSNE(n_components=2, random_state=0)
tsne_results = tsne.fit_transform(df[training_feat_names].values)
df_tsne = df.copy()
df_tsne['x-tsne'] = tsne_results[:,0]
df_tsne['y-tsne'] = tsne_results[:,1]

plt.style.use('fivethirtyeight')
fig,ax = plt.subplots(1,2, figsize = (10,6))
for cre_line in all_cre_lines:
    ax[0].scatter(df_tsne.loc[df_tsne['Cre_line'] == cre_line,  'x-tsne'], 
      df_tsne.loc[df_tsne['Cre_line'] == cre_line,  'y-tsne'],
      alpha = 0.5, s = 100, 
      label = cre_line) 
ax[0].set_xlabel('x-tsne')
ax[0].set_ylabel('y-tsne')

pca = PCA(n_components=2)
pca_result = pca.fit_transform(df[training_feat_names].values)
df_pca = df.copy()
df_pca['pca-one'] = pca_result[:,0]
df_pca['pca-two'] = pca_result[:,1] 

for cre_line in all_cre_lines:
    ax[1].scatter(df_pca.loc[df_tsne['Cre_line'] == cre_line,  'pca-one'], 
      df_pca.loc[df_pca['Cre_line'] == cre_line,  'pca-two'],
      alpha = 0.5, s = 100, 
      label = cre_line) 
ax[1].set_xlabel('pca1 (%.1f%s)'%(pca.explained_variance_ratio_[0]*100,'%'))
ax[1].set_ylabel('pca2 (%.1f%s)'%(pca.explained_variance_ratio_[1]*100,'%'))
plt.show()


svclassifier = SVC(kernel='rbf')  
svclassifier.fit(df[training_feat_names], df['Cre_line']) 
Cre_line_pred = svclassifier.predict(df[training_feat_names]) 
from sklearn.metrics import classification_report, confusion_matrix  
confustion_matrix_svm = confusion_matrix( df['Cre_line'], Cre_line_pred) 
print(classification_report( df['Cre_line'], Cre_line_pred))  
labels = sorted(list(set(df['Cre_line'].tolist())))
df_conf_svm = pd.DataFrame(confustion_matrix_svm, labels,
                  labels)

#fig,ax = plt.subplots()
#plt.style.use('ggplot')
#conf_svm = ax.imshow(confustion_matrix_svm,cmap='coolwarm');
#ax.grid('off')
#fig.colorbar(conf_svm, ax = ax)
#plt.show()

fig = plt.figure()
plt.style.use('fivethirtyeight')
sn.set(font_scale=1.4)#for label size
sn.heatmap(df_conf_svm, annot=True,annot_kws={"size": 16})# font size
fig.suptitle('Confusion Matrix for SVM classifier', fontsize = 14)
plt.show()

rnd_forest_classifier = RandomForestClassifier()
rnd_forest_classifier.fit(df[training_feat_names], df['Cre_line']) 
Cre_line_pred_rf = rnd_forest_classifier.predict(df[training_feat_names]) 
confustion_matrix_rnd_forest = confusion_matrix( df['Cre_line'], Cre_line_pred_rf) 

print(classification_report( df['Cre_line'], Cre_line_pred_rf))  
df_conf_rf = pd.DataFrame(confustion_matrix_rnd_forest, labels,
                  labels)
feature_imp = pd.Series(rnd_forest_classifier.feature_importances_,index=training_feat_names).sort_values(ascending=False)


fig = plt.figure()
plt.style.use('fivethirtyeight')
sn.set(font_scale=1.4)#for label size
sn.heatmap(df_conf_rf, annot=True,annot_kws={"size": 16})# font size
fig.suptitle('Confusion Matrix for Random Forest classifier', fontsize = 14)
plt.show()

fig = plt.figure()
sn.barplot(x=feature_imp, y=feature_imp.index)
plt.xlabel('Feature Importance Score')
plt.ylabel('Features')
plt.title("Visualizing Important Features")
plt.legend()
ax = plt.gca()
ax.tick_params(axis = 'y', which = 'major', labelsize = 10)
plt.show()