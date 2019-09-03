from unittest import TestCase
import os
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_passive,\
                                filter_feat_proto_active,adjust_param_bounds
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

def convert_stim_feat_to_dict(stim_features,key_name):
    model_feat_dict = {}
    for model_stim_,model_stim_feat in stim_features.items():
        for feat_name,feat_mean in model_stim_feat[key_name].items():
            model_feat_dict[model_stim_+'.'+ feat_name]= feat_mean[0] 
    return model_feat_dict
    

def convert_protocols_to_dict(stim_features,key_name):
    model_feat_dict = {}
    for model_stim_,model_stim_feat in stim_features.items():
        for feat_name,feat_mean in model_stim_feat[key_name][0].items():
            if 'sweep_filenames':
                continue
            model_feat_dict[model_stim_+'.'+ feat_name]= feat_mean
    return model_feat_dict

class Test_Model(TestCase):
    
    def setUp(self):
        self.mouse_spiny_id = '483101699'
        self.mouse_spiny_path = os.path.join(test_data_path,'mouse_spiny')
        
        self.mouse_aspiny_id = '481001895'
        self.mouse_aspiny_path = os.path.join(test_data_path,'mouse_aspiny')
        
        self.maxDiff = None
   
    # Check Stage0 stuff        
    def test_Stage0_parameters(self):
        
        # Mouse spiny
        cell_id_spiny = self.mouse_spiny_id 
        
        
        # Create the parameter bounds for the optimization
        model_params_handler_spiny = AllActive_Model_Parameters(cell_id_spiny)
        param_bounds_file = 'param_bounds_stage0.json'
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        model_params_spiny,model_params_release = model_params_handler_spiny.get_opt_params(param_bounds_path)
        model_mechs_spiny,model_mechs_release = model_params_handler_spiny.get_opt_mechanism(model_params_spiny,\
                            model_params_release,param_bounds_path)
        model_params_dict_spiny = convert_model_params_to_dict(model_params_spiny)
        
        mouse_spiny_stage0_params = os.path.join(self.mouse_spiny_path,'Stage0_parameters.json')
        mouse_spiny_stage0_mechs = os.path.join(self.mouse_spiny_path,'Stage0_mechanism.json')
        model_params_spiny_true = utility.load_json(mouse_spiny_stage0_params)
        model_mechs_spiny_true = utility.load_json(mouse_spiny_stage0_mechs)
        model_params_spiny_true_dict = convert_model_params_to_dict(model_params_spiny_true)
        
        self.assertEqual(model_params_spiny_true_dict,model_params_dict_spiny)
        self.assertEqual(model_mechs_spiny_true,model_mechs_spiny)
        
        # Mouse aspiny
        cell_id_aspiny = self.mouse_aspiny_id 
        
        # Create the parameter bounds for the optimization
        model_params_handler_aspiny = AllActive_Model_Parameters(cell_id_aspiny,\
                                                 swc_search_pattern = 'reconstruction.swc')
        param_bounds_file = 'param_bounds_stage0.json'
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        model_params_aspiny,model_params_release = model_params_handler_aspiny.\
                                            get_opt_params(param_bounds_path)
        model_mechs_aspiny,model_mechs_release = model_params_handler_aspiny.\
                                          get_opt_mechanism(model_params_aspiny,\
                                          model_params_release,param_bounds_path)
        model_params_dict_aspiny = convert_model_params_to_dict(model_params_aspiny)
        
        
        mouse_aspiny_stage0_params = os.path.join(self.mouse_aspiny_path,'Stage0_parameters.json')
        mouse_aspiny_stage0_mechs = os.path.join(self.mouse_aspiny_path,'Stage0_mechanism.json')
        model_params_aspiny_true = utility.load_json(mouse_aspiny_stage0_params)
        model_mechs_aspiny_true = utility.load_json(mouse_aspiny_stage0_mechs)
        model_params_aspiny_true_dict = convert_model_params_to_dict(model_params_aspiny_true)
        
        self.assertEqual(model_params_aspiny_true_dict,model_params_dict_aspiny)
        self.assertEqual(model_mechs_aspiny_true,model_mechs_aspiny)
    
    
    def test_Stage0_features(self):
        
        # Mouse spiny
        cell_id_spiny = self.mouse_spiny_id 
        
        nwb_handler_spiny = NWB_Extractor(cell_id_spiny,nwb_search_pattern=cell_id_spiny)
        acceptable_stimtypes = ['Long Square']
        ephys_dir = os.path.join(self.mouse_spiny_path,'mouse_spiny_ephys')
        ephys_data_path,stimmap_filename = \
                    nwb_handler_spiny.save_cell_data(acceptable_stimtypes,
                    ephys_dir=ephys_dir)
        feature_path = utility.locate_template_file(os.path.join('parameters',\
                            'feature_set_stage0.json'))
        train_features_spiny,_,_,train_protocols_spiny,_ = \
            nwb_handler_spiny.get_ephys_features(feature_path,ephys_data_path,
                                           stimmap_filename,filter_feat_proto_passive)
       
        train_features_dict_spiny = convert_stim_feat_to_dict(train_features_spiny,'soma')
        train_protocols_dict_spiny = convert_protocols_to_dict(train_protocols_spiny,'stimuli')
        mouse_spiny_stage0_features = os.path.join(self.mouse_spiny_path,'Stage0_features.json')
        mouse_spiny_stage0_protocols = os.path.join(self.mouse_spiny_path,'Stage0_protocols.json')
        model_features_spiny_true = utility.load_json(mouse_spiny_stage0_features)
        model_protocols_spiny_true = utility.load_json(mouse_spiny_stage0_protocols)
        model_features_spiny_true_dict = convert_stim_feat_to_dict(model_features_spiny_true,'soma')
        model_protocols_spiny_true_dict = convert_protocols_to_dict(model_protocols_spiny_true,
                                                              'stimuli')

        self.assertAlmostEqual(model_features_spiny_true_dict,train_features_dict_spiny)
        self.assertAlmostEqual(model_protocols_spiny_true_dict,train_protocols_dict_spiny)
        
        
        # Mouse aspiny
        cell_id_aspiny = self.mouse_aspiny_id 
        
        nwb_handler_aspiny = NWB_Extractor(cell_id_aspiny,nwb_search_pattern=cell_id_aspiny)
        acceptable_stimtypes = ['Long Square']
        ephys_dir = os.path.join(self.mouse_aspiny_path,'mouse_aspiny_ephys')
        ephys_data_path,stimmap_filename = \
                    nwb_handler_aspiny.save_cell_data(acceptable_stimtypes,
                    ephys_dir=ephys_dir)
        feature_path = utility.locate_template_file(os.path.join('parameters',\
                            'feature_set_stage0.json'))
        train_features_aspiny,_,_,train_protocols_aspiny,_ = \
            nwb_handler_aspiny.get_ephys_features(feature_path,ephys_data_path,
                                           stimmap_filename,filter_feat_proto_passive)
       
        train_features_dict_aspiny = convert_stim_feat_to_dict(train_features_aspiny,'soma')
        train_protocols_dict_aspiny = convert_protocols_to_dict(train_protocols_aspiny,'stimuli')
        mouse_aspiny_stage0_features = os.path.join(self.mouse_aspiny_path,'Stage0_features.json')
        mouse_aspiny_stage0_protocols = os.path.join(self.mouse_aspiny_path,'Stage0_protocols.json')
        model_features_aspiny_true = utility.load_json(mouse_aspiny_stage0_features)
        model_protocols_aspiny_true = utility.load_json(mouse_aspiny_stage0_protocols)
        model_features_aspiny_true_dict = convert_stim_feat_to_dict(model_features_aspiny_true,'soma')
        model_protocols_aspiny_true_dict = convert_protocols_to_dict(model_protocols_aspiny_true,
                                                              'stimuli')

        self.assertAlmostEqual(model_features_aspiny_true_dict,train_features_dict_aspiny)
        self.assertAlmostEqual(model_protocols_aspiny_true_dict,train_protocols_dict_aspiny)
        

    # Check Stage1 stuff        
    def test_Stage1_parameters(self):
        
        # Mouse spiny
        cell_id_spiny = self.mouse_spiny_id 
        
        # Create the parameter bounds for the optimization
        prev_params_file = 'Stage0_fit_spiny.json'
        model_params_handler_spiny = AllActive_Model_Parameters(cell_id_spiny,
                                                    prev_model_pattern=prev_params_file)
        param_bounds_file = 'param_bounds_stage1.json'
        
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        param_rule_func = adjust_param_bounds    
        model_params_spiny,model_params_release = model_params_handler_spiny.get_opt_params\
                                        (param_bounds_path,param_rule_func)
        model_mechs_spiny,model_mechs_release = model_params_handler_spiny.\
                            get_opt_mechanism(model_params_spiny,model_params_release,param_bounds_path)
        model_params_dict_spiny = convert_model_params_to_dict(model_params_spiny)
        mouse_spiny_stage1_params = os.path.join(self.mouse_spiny_path,'Stage1_parameters.json')
        mouse_spiny_stage1_mechs = os.path.join(self.mouse_spiny_path,'Stage1_mechanism.json')
        model_params_spiny_true = utility.load_json(mouse_spiny_stage1_params)
        model_mechs_spiny_true = utility.load_json(mouse_spiny_stage1_mechs)
        model_params_spiny_true_dict = convert_model_params_to_dict(model_params_spiny_true)
        
        self.assertEqual(model_params_spiny_true_dict,model_params_dict_spiny)
        self.assertEqual(model_mechs_spiny_true,model_mechs_spiny)
        
        # Mouse aspiny
        cell_id_aspiny = self.mouse_aspiny_id 
        
        
        # Create the parameter bounds for the optimization
        prev_params_file = 'Stage0_fit_aspiny.json'
        model_params_handler_aspiny = AllActive_Model_Parameters(cell_id_aspiny,\
                 swc_search_pattern = 'reconstruction.swc', prev_model_pattern=prev_params_file)
        param_bounds_file = 'param_bounds_stage1.json'
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        param_rule_func = adjust_param_bounds
        model_params_aspiny,model_params_release = model_params_handler_aspiny.\
                                    get_opt_params(param_bounds_path,param_rule_func)
        model_mechs_aspiny,model_mechs_release = model_params_handler_aspiny.\
                                          get_opt_mechanism(model_params_aspiny,\
                                          model_params_release,param_bounds_path)
        model_params_dict_aspiny = convert_model_params_to_dict(model_params_aspiny)
        
        model_params_dict_aspiny = convert_model_params_to_dict(model_params_aspiny)
        mouse_aspiny_stage1_params = os.path.join(self.mouse_aspiny_path,'Stage1_parameters.json')
        mouse_aspiny_stage1_mechs = os.path.join(self.mouse_aspiny_path,'Stage1_mechanism.json')
        model_params_aspiny_true = utility.load_json(mouse_aspiny_stage1_params)
        model_mechs_aspiny_true = utility.load_json(mouse_aspiny_stage1_mechs)
        model_params_aspiny_true_dict = convert_model_params_to_dict(model_params_aspiny_true)
        
        self.assertEqual(model_params_aspiny_true_dict,model_params_dict_aspiny)
        self.assertEqual(model_mechs_aspiny_true,model_mechs_aspiny)

    
    def test_Stage1_features(self):
        pass
    
    
    # Check Stage2 stuff        
    def test_Stage2_parameters(self):
        
        # Mouse spiny
        cell_id_spiny = self.mouse_spiny_id 
        
        # Create the parameter bounds for the optimization
        prev_params_file = 'Stage1_fit_spiny.json'
        model_params_handler_spiny = AllActive_Model_Parameters(cell_id_spiny,
                                                                prev_model_pattern=prev_params_file)
        param_bounds_file = 'param_bounds_stage2.json'
        
        param_bounds_path = utility.locate_template_file(os.path.join('parameters',\
                                                param_bounds_file))
        param_rule_func = adjust_param_bounds    
        model_params_spiny,model_params_release = model_params_handler_spiny.get_opt_params\
                                        (param_bounds_path,param_rule_func)
        model_mechs_spiny,model_mechs_release = model_params_handler_spiny.get_opt_mechanism(model_params_spiny,\
                            model_params_release,param_bounds_path)
        model_params_dict_spiny = convert_model_params_to_dict(model_params_spiny)
        mouse_spiny_stage2_params = os.path.join(self.mouse_spiny_path,'Stage2_parameters.json')
        mouse_spiny_stage2_mechs = os.path.join(self.mouse_spiny_path,'Stage2_mechanism.json')
        model_params_spiny_true = utility.load_json(mouse_spiny_stage2_params)
        model_mechs_spiny_true = utility.load_json(mouse_spiny_stage2_mechs)
        model_params_spiny_true_dict = convert_model_params_to_dict(model_params_spiny_true)
        
        self.assertEqual(model_params_spiny_true_dict,model_params_dict_spiny)
        
        for mech_key,mech_val in model_mechs_spiny.items():
            self.assertEqual(set(model_mechs_spiny_true[mech_key]),set(mech_val))

        # Mouse aspiny
        cell_id_aspiny = self.mouse_aspiny_id 
        
        
        # Create the parameter bounds for the optimization
        prev_params_file = 'Stage1_fit_aspiny.json'
        model_params_handler_aspiny = AllActive_Model_Parameters(cell_id_aspiny,\
                 swc_search_pattern = 'reconstruction.swc', prev_model_pattern=prev_params_file)
        param_bounds_file = 'param_bounds_stage2_mouse_aspiny.json'
        param_bounds_path = os.path.join(self.mouse_aspiny_path,param_bounds_file)
        param_rule_func = adjust_param_bounds
        model_params_aspiny,model_params_release = model_params_handler_aspiny.\
                                    get_opt_params(param_bounds_path,param_rule_func)
        model_mechs_aspiny,model_mechs_release = model_params_handler_aspiny.\
                                          get_opt_mechanism(model_params_aspiny,\
                                          model_params_release,param_bounds_path)
        model_params_dict_aspiny = convert_model_params_to_dict(model_params_aspiny)
        
        model_params_dict_aspiny = convert_model_params_to_dict(model_params_aspiny)
        mouse_aspiny_stage2_params = os.path.join(self.mouse_aspiny_path,'Stage2_parameters.json')
        mouse_aspiny_stage2_mechs = os.path.join(self.mouse_aspiny_path,'Stage2_mechanism.json')
        model_params_aspiny_true = utility.load_json(mouse_aspiny_stage2_params)
        model_mechs_aspiny_true = utility.load_json(mouse_aspiny_stage2_mechs)
        model_params_aspiny_true_dict = convert_model_params_to_dict(model_params_aspiny_true)
        
        self.assertEqual(model_params_aspiny_true_dict,model_params_dict_aspiny)
        for mech_key,mech_val in model_mechs_aspiny.items():
            self.assertEqual(set(model_mechs_aspiny_true[mech_key]),set(mech_val))
        
    def test_Stage2_features(self):
        pass

