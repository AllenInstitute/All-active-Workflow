import logging
from allensdk.core.cell_types_cache import CellTypesCache
import pandas as pd
from ateamopt.utils import utility
from ateamopt.analysis.optim_analysis import Optim_Analyzer
from functools import partial
from multiprocessing import Pool
import multiprocessing
import numpy as np
import os
from collections import defaultdict
from sklearn.preprocessing import StandardScaler   
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn import preprocessing
from sklearn.metrics import classification_report,\
                 confusion_matrix,accuracy_score  
from sklearn.model_selection import train_test_split  
from sklearn.utils.multiclass import unique_labels     
import seaborn as sns
import matplotlib.pyplot as plt
 
logger = logging.getLogger(__name__)

class Allactive_Classification(object):
    def __init__(self, param_file_list=None,metadata_file_list=None,
                     model_perf_filelist=None,me_cluster_data=None,
                     sdk_datapath=None,
                     ephys_file_list = None,
                     morph_file_list = None,
                     species=None):
        self.species = species
        self.param_file_list = param_file_list
        self.metadata_file_list = metadata_file_list
        self.sdk_datapath = sdk_datapath
        self.me_cluster_data = me_cluster_data
        self.model_perf_filelist = model_perf_filelist
        self.ephys_file_list = ephys_file_list
        self.morph_file_list = morph_file_list
    
    def allactive_param_data(self,repeat_params):
        analysis_handler = Optim_Analyzer()
        func = partial(analysis_handler.convert_aibs_param_to_dict,
                       repeat_params = repeat_params)
        p = Pool(multiprocessing.cpu_count())
        param_dict_list = p.map(func,self.param_file_list)
        p.close()
        p.join()
        param_df = pd.DataFrame(param_dict_list)
        return param_df
    
    def allactive_metadata(self,save_data=False):
        metadata_file_list = self.metadata_file_list
        metadata_list = []
        for metadata_file_ in metadata_file_list:
            metadata_list.append(utility.load_json(metadata_file_))
            
        metadata_df = pd.DataFrame(metadata_list)
        if save_data:
            self.save_class_data(metadata_df,'allactive_metadata.csv',
                        'allactive_metadatatype.csv')
        return metadata_df
        
    def morph_data(self,save_data=False):
        morph_file_list = self.morph_file_list
        morph_data_list = []
        for morph_file_ in morph_file_list:
            morph_dict = utility.load_json(morph_file_)
            cell_id = morph_file_.split('/')[-2]
            morph_dict['Cell_id'] = cell_id
            morph_data_list.append(morph_dict.copy())
         
        morph_df = pd.DataFrame(morph_data_list)   
        if save_data:
            self.save_class_data(morph_df,'morph_data.csv',
                        'morph_datatype.csv')
        return morph_df
    
    def sdk_data(self,save_data=False):
        if self.sdk_datapath and not os.path.exists(self.sdk_datapath):
            ctc = CellTypesCache()
            cells_allensdk = ctc.get_cells(species = [self.species],simple = False)
            cells_sdk_df = pd.DataFrame(cells_allensdk)
            cells_sdk_df.rename(columns={'specimen__id':'Cell_id'},inplace=True)
            cells_sdk_df['Cell_id'] = cells_sdk_df['Cell_id'].astype(str)
            if save_data:
                self.save_class_data(cells_sdk_df,self.sdk_datapath,
                                'sdk_datatype.csv')
        else:
            cells_sdk_df = self.read_class_data(self.sdk_datapath,
                                               'sdk_datatype.csv')
        return cells_sdk_df
    
    def me_data(self):
        me_cluster_data = self.me_cluster_data
        me_cluster_df = pd.read_csv(me_cluster_data,index_col=None)
        me_cluster_df.rename(columns={'specimen_id':'Cell_id'},inplace=True)
        me_cluster_df['Cell_id'] = me_cluster_df['Cell_id'].astype(str)
        return me_cluster_df
    
    def model_performance_data(self,save_data=False):
        metric_list_df = []
        model_perf_filelist = self.model_perf_filelist
        for metric_file_ in model_perf_filelist:
            cell_id= metric_file_.split('/')[-2]
            metric_list = utility.load_pickle(metric_file_)
            metric_dict = [{'hof_index' : ii, 'Feature_Avg_Train' : metric_list_['Feature_Average'],
                           'Feature_Avg_Generalization' : metric_list_['Feature_Average_Generalization'],
                           'Explained_Variance' : metric_list_['Explained_Variance'],
                           'Seed_Index' : metric_list_['Seed'],
                           'Cell_id' : cell_id} \
                           for ii,metric_list_ in enumerate(metric_list)]
            metric_list_df.extend(metric_dict)
        perf_metric_df = pd.DataFrame(metric_list_df)
        perf_metric_df['Explained_Variance'] *= 100
        
        if save_data:
            self.save_class_data(perf_metric_df,'performance_metric.csv',
                        'performance_datatype.csv')
        return perf_metric_df
    
    def ephys_data(self,save_data=False):
        select_features = ['AHP_depth',
                       'AP_amplitude_from_voltagebase',
                       'AP_width']
        ephys_file_list=self.ephys_file_list
        features_list = []
        for ephys_file_ in ephys_file_list:
            cell_id = ephys_file_.split('/')[-2]
            ephys_features_dict = {'Cell_id':cell_id}
            feature_dict = utility.load_json(ephys_file_)
            eFEL_features = self.get_eFEL_features(feature_dict,
                                               select_features)
            ephys_features_dict.update(eFEL_features)
            features_list.append(ephys_features_dict.copy())
        ephys_features_df = pd.DataFrame(features_list)
        
        return ephys_features_df
        
    @staticmethod
    def get_data_fields(data_path):
        if isinstance(data_path,pd.DataFrame):
            return list(data_path)
        else:
            if data_path.endswith('.json'):
                json_data = utility.load_json(data_path) 
                return list(json_data.keys())
            elif data_path.endswith('.csv'):
                csv_data = pd.read_csv(data_path, index_col = None)
                return list(csv_data)
        
        logger.debug('Not in .json,.csv or pandas dataframe format')
        return None
    
    @staticmethod
    def save_class_data(data,data_filename,datatype_filename):
        data.dtypes.to_frame('types').to_csv(datatype_filename)
        data.to_csv(data_filename, index=None)
    
    @staticmethod    
    def read_class_data(data_filename,datatype_filename):
        datatypes = pd.read_csv(datatype_filename)['types']
        data = pd.read_csv(data_filename, dtype=datatypes.to_dict())
        return data
        
    # Get cell ids for optimization from ME data
    @staticmethod        
    def get_cellid_for_opt(df_L,df_S,field='me_type',
                           non_std_morph_path=None,select_cre=None):
        df_L_select = df_L.loc[df_L[field].notnull(),['Cell_id',
                               'cre','dendrite_type',field]]
        if select_cre:
            df_L_select = df_L_select.loc[df_L_select['cre'].isin(select_cre),]
        
        df_L_cellids = df_L_select['Cell_id'].unique()
        df_S_cellids = df_S['Cell_id'].unique()
        
        if non_std_morph_path:
            non_std_morph_cells = pd.read_csv(non_std_morph_path,index_col=None)
            non_std_cell_ids = non_std_morph_cells.values
            df_S_cellids = np.append(df_S_cellids,non_std_cell_ids)
        cell_id_target = [cell_ for cell_ in df_L_cellids \
                          if cell_ not in df_S_cellids]
        df_L_select = df_L_select.loc[df_L_select['Cell_id'].isin(cell_id_target),]

        return df_L_select
        
    @staticmethod    
    def get_data_stat(data,field='Cre_line',agg_field='me_type'):
        data_grouped = data.groupby(field)
        grouped_stat = data_grouped[agg_field].agg(np.size)
        return grouped_stat
        
    @staticmethod    
    def get_eFEL_features(feature_dict,select_features):
        eFEL_features= defaultdict(list)
        for feat_stim,feat_val in feature_dict.items():
            eFEL_feat = feat_val['soma']
            for key_,val_ in eFEL_feat.items():
                if key_ in select_features:
                    eFEL_features[key_].append(val_[0])
                    
        eFEL_features = {feat_key:np.mean(feat_val) \
                 for feat_key,feat_val in eFEL_features.items()}
        return eFEL_features
    
    
    @staticmethod
    def prepare_data_clf(data,feature_fields,target_field,
                         least_pop=5):
        data_section = data.loc[:,feature_fields+[target_field]] 
        
        # drop any cell with target field nan
        data_section = data_section.dropna(axis=0, how = 'any',
                                           subset=[target_field])
        
        # filtering based on least populated class
        agg_data = data_section.groupby(target_field)[feature_fields[0]].\
                        agg(np.size).to_dict()
        filtered_targets = [key for key,val in agg_data.items() \
                            if val > least_pop]
        data_section = data_section.loc[data_section[target_field].\
                        isin(filtered_targets),]
        
        # drop any feature which is nan for any cells
        data_section = data_section.dropna(axis =1, how = 'any')
        revised_features = [feature_field for feature_field in \
                    list(data_section) if feature_field in feature_fields]
        X_df = data_section.loc[:,revised_features]
        y_df = data_section.loc[:,[target_field]]
        return X_df,y_df,revised_features
    
    @staticmethod    
    def SVM_classifier(X_df,y_df,feature_fields,
                       target_field,
                       plot_confusion_mat=False,
                       conf_mat_figname=None):
        
        
        clf = SVC(kernel='rbf')  
        svm_pipeline =  Pipeline([('scaler', StandardScaler()),
                                      ('svc', clf)])
        le = preprocessing.LabelEncoder()  
        y_df['label_encoder']= le.fit_transform(y_df[target_field])
        
        X_data = X_df.loc[:,feature_fields].values
        y_data = y_df['label_encoder'].values
        
        X_train, X_test, y_train, y_test = train_test_split(X_data, y_data,\
                        test_size=0.2, 
                        stratify=y_data, random_state=0)
        
        svm_pipeline.fit(X_train, y_train)
        y_pred_test = svm_pipeline.predict(X_test)
        confusion_matrix_svm = confusion_matrix(y_test, y_pred_test)
        svm_report = classification_report(y_test, y_pred_test)
        score = accuracy_score(y_test, y_pred_test)*100
        
        if plot_confusion_mat:
            classes = le.inverse_transform(unique_labels(y_test, \
                                            y_pred_test))
            df_conf_svm = pd.DataFrame(confusion_matrix_svm, classes,
                  classes)
            fig = plt.figure()
            sns.set(style="darkgrid", font_scale=1)
            ax = sns.heatmap(df_conf_svm, annot=True, fmt="d")# font size
            fig.suptitle('Accuracy: %s %%'%score, fontsize = 14)
            ax.set_ylabel('True Label')
            ax.set_xlabel('Predicted Label')
            fig.savefig(conf_mat_figname,bbox_inches='tight')
            
        return score,confusion_matrix_svm,svm_report
            
    
    