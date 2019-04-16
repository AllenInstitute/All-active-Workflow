from unittest import TestCase
import os
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_passive
#import nose.tools as nt


test_data_path = os.path.join(os.path.dirname(
                    os.path.abspath(__file__)),'test_data')
        

class Test_Model(TestCase):
    def test_hello(self):
        s = 1
        self.assertTrue(s,1)

    # Check Stage0 stuff        
    def test_Stage0_parameters(self):
        
        # Mouse spiny
        cell_id = '483101699'
        mouse_spiny_path = os.path.join(test_data_path,'mouse_spiny')
        
        # Create the parameter bounds for the optimization
        model_params_handler = AllActive_Model_Parameters(cell_id)
        param_bounds_file = 'param_bounds_stage0.json'
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        model_params,model_params_release = model_params_handler.get_opt_params(param_bounds_path)
        model_mechs,model_mechs_release = model_params_handler.get_opt_mechanism(model_params,\
                            model_params_release,param_bounds_path)
        
        mouse_spiny_stage0_params = os.path.join(mouse_spiny_path,'Stage0_parameters.json')
        mouse_spiny_stage0_mechs = os.path.join(mouse_spiny_path,'Stage0_mechanism.json')
        model_params_true = utility.load_json(mouse_spiny_stage0_params)
        model_mechs_true = utility.load_json(mouse_spiny_stage0_mechs)
        
        self.assertListEqual(model_params_true,model_params)
        self.assertListEqual(model_mechs_true,model_mechs)
    
    
def test_Stage0_features():
    
    # Mouse spiny
#    cell_id = '483101699'
#    
#    nwb_handler = NWB_Extractor(cell_id)
#    acceptable_stimtypes = ['Long Square']
#    ephys_data_path,stimmap_filename = \
#                nwb_handler.save_cell_data(acceptable_stimtypes)
#    feature_path = utility.locate_template_file(os.path.join('parameters',\
#                        'feature_set_stage0.json'))
#    filter_rule_func = filter_feat_proto_passive
#    features_write_path,untrained_features_write_path,all_features_write_path,\
#        protocols_write_path,all_protocols_write_path = \
#        nwb_handler.get_ephys_features(feature_path,ephys_data_path,
#                                       stimmap_filename,filter_rule_func)
    pass

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
