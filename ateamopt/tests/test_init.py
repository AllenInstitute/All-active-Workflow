from unittest import TestCase
import os
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_passive




class Test_Model(TestCase):
    def test_hello(self):
        s = 1
        self.assertTrue(s,1)

test_data_path = os.path.join(os.path.dirname(
                    os.path.abspath(__file__)),'test_data')
        
# Check Stage0 stuff        
def test_Stage0_parameters():
    
    # Mouse spiny
    cell_id = '483101699'
    
    # Create the parameter bounds for the optimization
    model_params_handler = AllActive_Model_Parameters(cell_id)
    param_bounds_file = 'param_bounds_stage0.json'
    param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                            param_bounds_file))
    model_params,model_params_release = model_params_handler.get_opt_params(param_bounds_path)


def test_Stage0_features():
    
    # Mouse spiny
    cell_id = '483101699'
    nwb_handler = NWB_Extractor(cell_id)
    acceptable_stimtypes = ['Long Square']
    ephys_data_path,stimmap_filename = \
                nwb_handler.save_cell_data(acceptable_stimtypes)
    feature_path = utility.locate_template_file(os.path.join('parameters',\
                        'feature_set_stage0.json'))
    filter_rule_func = filter_feat_proto_passive
    features_write_path,untrained_features_write_path,all_features_write_path,\
        protocols_write_path,all_protocols_write_path = \
        nwb_handler.get_ephys_features(feature_path,ephys_data_path,
                                       stimmap_filename,filter_rule_func)

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
