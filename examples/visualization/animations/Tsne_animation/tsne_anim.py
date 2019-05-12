from ateamopt.utils import utility
import pandas as pd
import numpy as np

import os
from ateamopt.animation import animation_module
from ateamopt.analysis.allactive_classification \
                import Allactive_Classification as aa_clf


def main():
    clf_handler = aa_clf()
    mouse_data_file = 'data/Mouse_class_data.csv'
    mouse_datatype_file = 'data/Mouse_class_datatype.csv'
    mouse_data_df = clf_handler.read_class_data(mouse_data_file,mouse_datatype_file)
    cre_data = mouse_data_df.loc[mouse_data_df.hof_index==0,\
                                    ['Cell_id','Cre_line']]
    
    all_active_param_file = 'data/allactive_params.csv'
    param_datatype_file = 'data/allactive_paramsdatatype.csv'
    hof_param_data = clf_handler.read_class_data(all_active_param_file,
                                                 param_datatype_file)
    
    max_hof_index = np.max(hof_param_data.hof_index.values)
    target_field = 'Cre_line'
    
    files = []
    lims = [-100,100]
    for hof_ind in range(max_hof_index+1):
        param_data = hof_param_data.loc[hof_param_data.hof_index<= \
                                        hof_ind,]
        param_data = param_data.drop(labels='hof_index',axis=1)
        df_tsne_p_cre = pd.merge(param_data,cre_data,how='left',
                            on='Cell_id')
    
        p_X_df,p_y_df,revised_features = clf_handler.prepare_data_clf\
                            (df_tsne_p_cre,list(param_data),target_field,
                             least_pop=5*(hof_ind+1))
        tsne_features = [feature_ for feature_ in revised_features \
                     if feature_ != 'Cell_id']                  
        anim_prefix = 'TSNE_figures'
        fname = 'tsne_{}.jpeg'.format(hof_ind)
        file_path = os.path.join(anim_prefix,fname)
        utility.create_filepath(file_path)
        clf_handler.tsne_embedding(p_X_df,p_y_df,tsne_features,
               target_field,marker_field=target_field,
               figname=file_path,
           figtitle = 'TSNE : Model parameters-{} models per cell'.format(hof_ind+1),
           xlims=lims,ylims=lims)
        files.append(file_path)
        
    anim_handler = animation_module.Animation(movie_name = \
                                      'Tsne_anim.gif')
    anim_handler.make_gif(files)    
        
    
if __name__ == '__main__':
    main()
