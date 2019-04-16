from unittest import TestCase
import os
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_passive
#import nose.tools as nt


test_data_path = os.path.join(os.path.dirname(
                    os.path.abspath(__file__)),'test_data')
        
def convert_model_params_to_dict(model_params_opt):
    model_param_dict = {}
    for model_param_ in model_params_opt:
        if model_param_.get('bounds'):
            model_param_dict[model_param_['param_name']+'.'+ \
                 model_param_['sectionlist']]= model_param_['bounds'] 
             
        elif model_param_.get('value'):
            model_param_dict[model_param_['param_name']+'.'+ \
                             model_param_.get('sectionlist','all')]= model_param_['value'] 
    
    return model_param_dict

def convert_stim_feat_to_dict(stim_features):
    model_param_dict = {}
    for model_param_ in model_params_opt:
        if model_param_.get('bounds'):
            model_param_dict[model_param_['param_name']+'.'+ \
                 model_param_['sectionlist']]= model_param_['bounds'] 
             
        elif model_param_.get('value'):
            model_param_dict[model_param_['param_name']+'.'+ \
                             model_param_.get('sectionlist','all')]= model_param_['value'] 
    
    return model_param_dict
    
class Test_Model(TestCase):
    
    def setUp(self):
        self.mouse_spiny_id = '483101699'
        self.mouse_spiny_path = os.path.join(test_data_path,'mouse_spiny')
   
    # Check Stage0 stuff        
    def test_Stage0_parameters(self):
        
        # Mouse spiny
        cell_id = self.mouse_spiny_id 
        
        
        # Create the parameter bounds for the optimization
        model_params_handler = AllActive_Model_Parameters(cell_id)
        param_bounds_file = 'param_bounds_stage0.json'
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        model_params,model_params_release = model_params_handler.get_opt_params(param_bounds_path)
        model_mechs,model_mechs_release = model_params_handler.get_opt_mechanism(model_params,\
                            model_params_release,param_bounds_path)
        model_params_dict = convert_model_params_to_dict(model_params)
        
        
        mouse_spiny_stage0_params = os.path.join(self.mouse_spiny_path,'Stage0_parameters.json')
        mouse_spiny_stage0_mechs = os.path.join(self.mouse_spiny_path,'Stage0_mechanism.json')
        model_params_true = utility.load_json(mouse_spiny_stage0_params)
        model_mechs_true = utility.load_json(mouse_spiny_stage0_mechs)
        model_params_true_dict = convert_model_params_to_dict(model_params_true)
        
        self.assertEqual(model_params_true_dict,model_params_dict)
        self.assertEqual(model_mechs_true,model_mechs)
    
    
    def test_Stage0_features(self):
        
        # Mouse spiny
        cell_id = self.mouse_spiny_id 
        
        nwb_handler = NWB_Extractor(cell_id,nwb_search_pattern=cell_id)
        acceptable_stimtypes = ['Long Square']
        ephys_dir = os.path.join(self.mouse_spiny_path,'mouse_spiny_ephys')
        ephys_data_path,stimmap_filename = \
                    nwb_handler.save_cell_data(acceptable_stimtypes,
                    ephys_dir=ephys_dir)
        feature_path = utility.locate_template_file(os.path.join('parameters',\
                            'feature_set_stage0.json'))
        train_features,_,_,train_protocols,_ = \
            nwb_handler.get_ephys_features(feature_path,ephys_data_path,
                                           stimmap_filename,filter_feat_proto_passive)
        mouse_spiny_stage0_features = os.path.join(self.mouse_spiny_path,'Stage0_features.json')
        mouse_spiny_stage0_protocols = os.path.join(self.mouse_spiny_path,'Stage0_protocols.json')
        model_params_true = utility.load_json(mouse_spiny_stage0_features)
        model_mechs_true = utility.load_json(mouse_spiny_stage0_protocols)
    

def test_Stage0_mechanisms():
    pass

def test_Stage0_protocols():
    pass
    


# Check Stage1 stuff        
def test_Stage1_parameters():
    pass

def test_Stage1_features():
    pass

def test_Stage1_mechanisms():
    pass

def test_Stage1_protocols():
    pass
    
# Check Stage2 stuff        
def test_Stage2_parameters():
    pass

def test_Stage2_features():
    pass

def test_Stage2_mechanisms():
    pass

def test_Stage2_protocols():
    pass
