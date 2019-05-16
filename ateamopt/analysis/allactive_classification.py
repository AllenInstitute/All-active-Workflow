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
from sklearn.manifold import TSNE
from sklearn.svm import SVC
from sklearn import preprocessing
from sklearn.metrics import classification_report,\
                 confusion_matrix,accuracy_score  
from sklearn.model_selection import train_test_split  
from sklearn.utils.multiclass import unique_labels     
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

 
logger = logging.getLogger(__name__)

marker_list = list(Line2D.markers.keys())[2:-4]

def reformat_legend(axes,sorted_class,color_list):
    
    all_handles,all_labels,all_colors = [], [], []
    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            handles,labels = axes[i][j].get_legend_handles_labels()
            for handle_,label_ in zip(handles,labels):
                label_split = label_.split('.')[0]
                if  label_split not in all_labels:
                    all_labels.append(label_split)
                    all_handles.append(handle_)
                    all_colors.append(color_list[sorted_class.index(label_split)])
    
    dummy_handles = []
    for handle_,label_,color_ in zip(all_handles,all_labels,all_colors):
        h = plt.scatter([],[],marker = 'o',color = color_,
                      label = label_)
        dummy_handles.append(h)
    
    return dummy_handles,all_labels


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
    
    def allactive_param_data(self,repeat_params,save_data=False):
        analysis_handler = Optim_Analyzer()
        func = partial(analysis_handler.convert_aibs_param_to_dict,
                       repeat_params = repeat_params)
        p = Pool(multiprocessing.cpu_count())
        param_dict_list = p.map(func,self.param_file_list)
        p.close()
        p.join()
        param_df = pd.DataFrame(param_dict_list)
        if save_data:
            self.save_class_data(param_df,'allactive_params.csv',
                        'allactive_paramsdatatype.csv')
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
    def get_cellid_for_opt(df_L,df_S,target_field='cre',
                           addl_target = 'me_type',
                           non_std_morph_path=None,select_types=None):
        info_fields = ['cre','dendrite_type']
        select_fields = list(set([target_field, addl_target]+info_fields))
        select_fields = [select_field_ for select_field_ in select_fields \
                         if select_field_ is not None]
        df_L = df_L.dropna(subset=select_fields,how='any',axis=0)
        df_L_select = df_L.loc[:,['Cell_id']+select_fields]
        if select_types:
            df_L_select = df_L_select.loc[df_L_select[target_field].\
                                          isin(select_types),]
        
        df_L_cellids = df_L_select['Cell_id'].unique()
        df_S_cellids = df_S['Cell_id'].unique()
        
        if non_std_morph_path:
            non_std_morph_cells = pd.read_csv(non_std_morph_path,index_col=None)
            non_std_cell_ids = non_std_morph_cells.values
            df_S_cellids = np.append(df_S_cellids,non_std_cell_ids)
        cell_id_target = [cell_ for cell_ in df_L_cellids \
                          if cell_ not in df_S_cellids]
        df_L_select = df_L_select.loc[df_L_select['Cell_id'].\
                                      isin(cell_id_target),]
        df_L_select = df_L_select.reset_index(drop=True)
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
                         property_fields = [],
                         least_pop=5):
        
        data = data.loc[:,~data.columns.duplicated()]
        data_section = data.loc[:,feature_fields+property_fields+\
                                [target_field]] 
        
        # drop any cell with target field nan
        data_section = data_section.dropna(axis=0, how = 'any',
                               subset=[target_field] + property_fields)
        
        # filtering based on least populated class
        agg_data = data_section.groupby(target_field)[feature_fields[0]].\
                        agg(np.size).to_dict()
        filtered_targets = [key for key,val in agg_data.items() \
                            if val > least_pop]
        data_section = data_section.loc[data_section[target_field].\
                        isin(filtered_targets),]
        
        # drop any feature which is nan for any cells
        data_section = data_section.dropna(axis =1,
                                           how = 'any')
        revised_features = [feature_field for feature_field in \
                    list(data_section) if feature_field in feature_fields]
        X_df = data_section.loc[:,revised_features + property_fields]
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
                        test_size=0.3, 
                        stratify=y_data, random_state=0)
        
        svm_pipeline.fit(X_train, y_train)
        y_pred_test = svm_pipeline.predict(X_test)
        confusion_matrix_svm = confusion_matrix(y_test, y_pred_test)
        svm_report = classification_report(y_test, y_pred_test)
        score = accuracy_score(y_test, y_pred_test)*100
        classes = le.inverse_transform(unique_labels(y_test, \
                                        y_pred_test))
        
        df_conf_svm = pd.DataFrame(confusion_matrix_svm, classes,
              classes)
        df_conf_svm=df_conf_svm.div(df_conf_svm.sum(axis=1),axis=0)
        
        if plot_confusion_mat:
            sns.set(style="darkgrid", font_scale=1)
            fig = plt.figure()
            cmap = sns.cubehelix_palette(light=1, as_cmap=True)
            ax = sns.heatmap(df_conf_svm,cmap=cmap,alpha=0.8)
            fig.suptitle('Accuracy: %s %%'%score, fontsize = 14)
            ax.set_ylabel('True Label')
            ax.set_xlabel('Predicted Label')
            fig.savefig(conf_mat_figname,bbox_inches='tight')
            plt.close(fig)
        return score,confusion_matrix_svm,svm_report,df_conf_svm
            
    
    @staticmethod
    def tsne_embedding(X_df,y_df,feature_fields,
                       target_field,marker_field='Cell_id',
                       col_var = None,figname='tsne_plot.pdf',
                       figtitle = '',
                       **kwargs):
        
        tsne_pipeline = Pipeline([('scaling', StandardScaler()), \
                     ('tsne',TSNE(n_components=2,random_state =0))])
            
        X_data = X_df.loc[:,feature_fields].values
        tsne_results = tsne_pipeline.fit_transform(X_data)
        
        data = pd.concat([X_df,y_df],axis=1)
        data['x-tsne'] = tsne_results[:,0]
        data['y-tsne'] = tsne_results[:,1]
        data['hue_marker'] = data.apply(lambda x : \
                    x[target_field] +'.'+ x[marker_field], axis =1)
        
        sorted_target = sorted(list(data[target_field].unique()))
        sorted_marker = sorted(list(data[marker_field].unique()))
        sorted_hue_marker = sorted(list(data['hue_marker'].unique()))
        sorted_col = sorted(list(data[col_var].unique())) \
                    if col_var else None
        
        current_palette = sns.color_palette()
        color_list = sns.color_palette(current_palette,len(sorted_target))
        color_df_list = [color_list[sorted_target.index(hue_.split('.')[0])] \
                         for hue_ in sorted_hue_marker]
        
        marker_df_list = [marker_list[sorted_marker.index(hue_.\
                           split('.')[1])%len(marker_list)] \
                         for hue_ in sorted_hue_marker]
        
        if target_field == marker_field:
            marker_df_list = ['o']*len(marker_df_list)
        
        sns.set(style="darkgrid", font_scale=1.5)
        g = sns.FacetGrid(data,
                          col=col_var,col_order=sorted_col, 
                          hue='hue_marker',hue_order = sorted_hue_marker,
                          palette=color_df_list, hue_kws=dict(marker=marker_df_list),                      
                          height=9,aspect= 1.2)
        g = (g.map(plt.scatter, "x-tsne", "y-tsne", alpha=.8,s=80))
         
        axes = g.axes
        dummy_handles,all_labels = reformat_legend(axes,sorted_target,color_list)
        g.fig.tight_layout(rect=[0, 0.1, 1, 0.95])
        g.fig.legend(handles = dummy_handles, labels = all_labels,
                     loc = 'lower center', ncol = 6)
        figtitle += ' (n=%s)'%len(data['Cell_id'].unique())
        g.fig.suptitle(figtitle)
        
#        if kwargs.get('xlims'):
#            for i in range(axes.shape[0]):
#                for j in range(axes.shape[1]):
#                    axes[i,j].set_xlim(kwargs['xlims'])
#        if kwargs.get('ylims'):
#             for i in range(axes.shape[0]):
#                for j in range(axes.shape[1]):
#                    axes[i,j].set_ylim(kwargs['ylims'])
        
        g.fig.savefig(figname, dpi= 80,bbox_inches='tight')
        plt.close(g.fig)
        
    @staticmethod
    def get_fi_slope(stim_fi, spike_fi):
        spiking_stim_idx = [ix for ix,rate in enumerate(spike_fi)\
                               if rate > 0]
        stim_fi = np.array(stim_fi)[spiking_stim_idx]
        spike_fi = np.array(spike_fi)[spiking_stim_idx]
        x = np.array(stim_fi)*1e3
        y = np.array(spike_fi)
        A = np.vstack([x, np.ones_like(x)]).T
        m, _ = np.linalg.lstsq(A, y)[0]
        return m

    @staticmethod
    def get_fi_intercept(stim_fi, spike_fi):
        nonspiking_stim_idx = [ix for ix,rate in enumerate(spike_fi) if rate ==0]
        non_spiking_stim = [stim for idx,stim in enumerate(stim_fi) \
                if idx in nonspiking_stim_idx and stim > 0]
        intercept = non_spiking_stim[-1]
        return intercept*1e3  
    
    
    def get_fI_prop(self,fi_data_file):
        if type(fi_data_file) is tuple:
            cell_id = fi_data_file[0].split('_')[-1].split('.')[0]
            data_dict = {'Cell_id' : cell_id}
            for fi_file in fi_data_file:
                fi_data = utility.load_pickle(fi_file)
                stim_key = [key for key in fi_data.keys() \
                                if 'stim' in key][0]
                freq_key = [key for key in fi_data.keys() \
                                if 'freq' in key][0]
                data_type = stim_key.split('_',1)[1]
                slope_ = self.get_fi_slope(fi_data[stim_key],fi_data[freq_key])
                icpt_ = self.get_fi_intercept(fi_data[stim_key],fi_data[freq_key])
                data_dict['slope_{}'.format(data_type)] = slope_
                data_dict['intercept_{}'.format(data_type)] = icpt_
        else:
            fi_data = utility.load_pickle(fi_data_file)
            cell_id = fi_data_file.split('_')[-1].split('.')[0]
            data_dict = {'Cell_id' : cell_id}
            stim_keys = sorted([key for key in fi_data.keys() \
                                if 'stim' in key])
            freq_keys = sorted([key for key in fi_data.keys() \
                                if 'freq' in key])
            for stim,freq in zip(stim_keys,freq_keys):
                data_type = stim.split('_',1)[1]
                slope_ = self.get_fi_slope(fi_data[stim],fi_data[freq])
                icpt_ = self.get_fi_intercept(fi_data[stim],fi_data[freq])
                data_dict['slope_{}'.format(data_type)] = slope_
                data_dict['intercept_{}'.format(data_type)] = icpt_
                
        return data_dict
        
        
        
        
        
        