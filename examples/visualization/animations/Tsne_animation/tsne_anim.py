from ateamopt.utils import utility
import pandas as pd
from sklearn.preprocessing import StandardScaler   
from sklearn.pipeline import Pipeline 
from sklearn.manifold import TSNE
import numpy as np
from matplotlib.lines import Line2D
import seaborn as sns
import matplotlib.pyplot as plt
import os
from ateamopt.animation import animation_module

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

def tsne_embedding(param_df,hue="Cre_line", 
            figname = 'tsne_%s.pdf'%np.random.randn(), 
            figtitle = '',drop_apical = False,
            force_hue = False,
            col_var = 'Layer',
            ignore_metadata_fields = None,
            include_metadata_fields = None):
    
    param_df = param_df.dropna(axis=1, how = 'all')
    
   # Get the conductance parameters
    if ignore_metadata_fields:
        param_df_var = param_df.loc[:,[col_name for col_name in \
                      list(param_df) if col_name not in ignore_metadata_fields]]
    elif include_metadata_fields:
        param_df_var = param_df.loc[:,[col_name for col_name in \
                      list(param_df) if col_name in include_metadata_fields]]
        param_df_var = param_df_var.dropna(axis=1, how = 'any')
    if drop_apical:
        cols = [c for c in param_df_var.columns if 'apic' not in c]
        param_df_var = param_df_var[cols]
    
    param_df_var = param_df_var.dropna(axis=0, how = 'any')
    
    
    X_cond = param_df_var.values
    tsne_pipeline = Pipeline([('scaling', StandardScaler()), \
                             ('tsne',TSNE(n_components=2,random_state =1))])
    
    tsne_results = tsne_pipeline.fit_transform(X_cond)
    param_df_tsne = param_df.loc[param_df_var.index]
    param_df_tsne[hue] = param_df_tsne[hue].astype(str)
    param_df_tsne['x-tsne'] = tsne_results[:,0]
    param_df_tsne['y-tsne'] = tsne_results[:,1]
    param_df_tsne['hue_id'] = param_df_tsne.apply(lambda x : \
                x[hue] +'.'+ x.Cell_id, axis =1)
    
    sorted_class = sorted(list(param_df_tsne[hue].unique()))
    sorted_cell_id = sorted(list(param_df_tsne.Cell_id.unique()))
    sorted_hue_id = sorted(list(param_df_tsne['hue_id'].unique()))
    if col_var:
        sorted_col = sorted(list(param_df_tsne[col_var].unique()))
    else:
        sorted_col = None
    current_palette = sns.color_palette()
    color_list = sns.color_palette(current_palette,len(sorted_class))
    color_df_list = [color_list[sorted_class.index(hue_.split('.')[0])] \
                     for hue_ in sorted_hue_id]
    marker_df_list = [marker_list[sorted_cell_id.index(hue_.\
                       split('.')[1])%len(marker_list)] \
                     for hue_ in sorted_hue_id]
    
    sns.set(style="darkgrid", font_scale=1)
    g = sns.FacetGrid(param_df_tsne,
                      col=col_var,col_order=sorted_col, 
                      hue='hue_id',hue_order = sorted_hue_id,
                      palette=color_df_list, hue_kws=dict(marker=marker_df_list),                      
                      height=9,aspect= 1.2)
    g = (g.map(plt.scatter, "x-tsne", "y-tsne", alpha=.8,s=80))
     
    axes = g.axes
    
    dummy_handles,all_labels = reformat_legend(axes,sorted_class,color_list)
    if not drop_apical:
        g.fig.legend(handles = dummy_handles, labels = all_labels,
                         loc = 'center right',fontsize=12)
        g.fig.tight_layout(rect=[0, 0.1, .82, 0.95])
    else:
        g.fig.legend(handles = dummy_handles, labels = all_labels,
                 loc = 'lower center', ncol = 6)
        g.fig.tight_layout(rect=[0, 0.1, 1, 0.95])
    
    figtitle += ' (# of cells =%s)'%len(sorted_cell_id)
    g.fig.suptitle(figtitle)
    g.fig.savefig(figname, dpi= 80)
    plt.close(g.fig)

def select_models_on_tolerance(df,tol_val):
    best_models_df = df.loc[df.hof_index == 0,:]
    cre_grouped = best_models_df.groupby('Cre_line')
    cre_grouped_metrics = cre_grouped['Feature_Avg_Train'].agg(np.mean)
    cre_grouped_metrics_dict = dict(cre_grouped_metrics)
    
    cre_grouped_metrics_tol = {key : (1+tol_val)*val for key,val in \
                               cre_grouped_metrics_dict.items()}
    df_tol = []
    df_hof = df.loc[df.hof_index == 0,] # Atleast have representative per cell
    for cre,fa_tol in cre_grouped_metrics_tol.items():
        df_temp = df.loc[(df.Cre_line == cre) & (df.Feature_Avg_Train <= fa_tol),]
        df_tol.append(df_temp)
        
    df_tol = pd.concat([pd.concat(df_tol),df_hof])
    df_tol = df_tol.drop_duplicates()
    return df_tol,cre_grouped_metrics_dict

def broad_cre_class(cre_var):
    if cre_var.startswith('Htr3a'):
        return 'Htr3a'
    elif cre_var.startswith('Pvalb'):
        return 'Pvalb'
    elif cre_var.startswith('Sst'):
        return 'Sst'
    else:
        return cre_var

#def main():
mouse_classification_data = 'Mouse_class_data.csv'
mouse_data_df = pd.read_csv(mouse_classification_data, 
                        dtype={'Cell_id': str},index_col = 0)
tol_level = .2
mouse_data_df,cre_grouped_metrics = \
            select_models_on_tolerance(mouse_data_df,tol_level)
#    mouse_data_df['Cre_line'] = mouse_data_df['Cre_line'].apply(broad_cre_class)
metadata_fields = utility.load_json('metadata_fields.json')
max_hof_index = 40
filter_feat = 'Dendrite_type'; filter_val = 'spiny'
drop_apical = filter_val == 'aspiny'
mouse_data_filtered = mouse_data_df.loc[mouse_data_df[filter_feat] == filter_val,]
mouse_data_hof = mouse_data_filtered

files = []
#for hof_ind in range(max_hof_index):
#    mouse_data_hof = mouse_data_filtered.loc[mouse_data_filtered['hof_index'] \
#                     <= hof_ind,]
#    anim_prefix = 'GA_figures_%s'%filter_val
#    fname = 'tsne_{}.jpeg'.format(hof_ind)
#    file_path = os.path.join(anim_prefix,fname)
#    utility.create_filepath(file_path)
tsne_embedding(mouse_data_hof,hue="Cre_line", 
        figname = file_path, 
        col_var = 'Species',
        figtitle = 'Mouse %s: Hall-of-fame index %s'%(filter_val,hof_ind),
        ignore_metadata_fields = metadata_fields,
        drop_apical = drop_apical)

#files.append(file_path)
#
#anim_handler = animation_module.Animation(movie_name = \
#                                  'Tsne_anim_{}.gif'.format(filter_val))
#anim_handler.make_gif(files)    
    
    
#if __name__ == '__main__':
#    main()
