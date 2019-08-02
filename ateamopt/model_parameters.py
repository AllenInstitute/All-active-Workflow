
from ateamopt.utils import utility
import re
import logging
from collections import defaultdict

logger = logging.getLogger(__name__)

class AllActive_Model_Parameters(object):

    def __init__(self,cell_id,temp=34,v_init=-80,swc_path=None,\
                 swc_search_pattern='cell_types',
                 released_aa_model_pattern='fit_aa_parameters',
                 prev_model_pattern='fit_opt',model_param_extension = '.json'):

        self.cell_id = cell_id
        self.v_init = v_init
        self.temperature = temp

        self._swc_path = swc_path
        if not swc_path and swc_search_pattern:
            swc_dir =  utility.get_filepath_for_exten(exten = '.swc')
            try:
                self._swc_path = [str_path for str_path in swc_dir if \
                                      swc_search_pattern in str_path][0]
            except:
                pass

        model_path_list = utility.get_filepath_for_exten(model_param_extension)

        # Parameter file for a released all-active model
        #(has to match name fit_aa_parameters.json)
        release_param_path_ = [str_path for str_path in model_path_list if \
                                   released_aa_model_pattern in str_path]
        try:
            self.release_param_path = release_param_path_[0]
        except:
            self.release_param_path = None

        # Parameter file from previous stage
        prev_stage_model_path_ = [str_path for str_path in model_path_list if \
                                   prev_model_pattern in str_path]
        try:
            self.prev_stage_model_path = prev_stage_model_path_[0]
        except:
            self.prev_stage_model_path = None


        self.no_apical = utility.check_swc_for_apical(self.swc_path) \
                if self.swc_path else False

    @property
    def swc_path(self):
        return self._swc_path


    @staticmethod
    def group_params(params_dict):
        active_params, Ih_params, passive_params = [],[],[]
        all_params_dict = {}

        ena_sect, ek_sect = [],[]

        for param_name,param_dict in params_dict.items():
            if re.search('Na', param_name, re.IGNORECASE):           # Use re.IGNORECASE if needed
                for sect in param_dict['section']:
                    ena_sect.append(sect)
            elif re.search('K', param_name, re.IGNORECASE):
                for sect in param_dict['section']:         # Use re.IGNORECASE if needed
                    ek_sect.append(sect)

        ena_sect = list(set(ena_sect))
        ek_sect = list(set(ek_sect))


        for param_name,param_item in params_dict.items():
            if 'mechanism' not in param_item.keys():
                passive_params.append(param_name)
            elif param_item['mechanism'] == 'Ih' or  param_item['mechanism'] == 'HCN':
                Ih_params.append(param_name)
            else:
                active_params.append(param_name)
            all_params_dict[param_name] = param_item

        return active_params, Ih_params, passive_params,all_params_dict,ena_sect,ek_sect


    def get_opt_params(self,param_bounds_path,adjust_param_rule = None):

        section_map = utility.bpopt_section_map
        params_dict = utility.load_json(param_bounds_path)
        _, _,_,all_params_dict,\
                ena_sect,ek_sect = self.group_params(params_dict)

        model_params_opt = list()

        # Create parameter file from the initialization bounds for optimization

        for param_name,param_dict in all_params_dict.items():
            for sect in param_dict['section']:

                 if self.no_apical and sect == 'apic': # if no apical dendrite in morphology
                     continue

                 iter_dict = {'param_name': param_name}
                 iter_dict['sectionlist'] = section_map[sect]
                 iter_dict['type'] = 'section'
                 iter_dict['dist_type'] = 'uniform'
                 try:
                     iter_dict['mech'] = param_dict['mechanism']
                     iter_dict['type'] = 'range'
                 except:
                     pass
                 iter_dict['bounds'] = param_dict['bounds'][sect]
                 model_params_opt.append(iter_dict)


        # Adjust parameter bounds from previous stage

        if self.prev_stage_model_path and adjust_param_rule:
            model_params_prev = self.load_params_prev_stage(section_map)
            for model_param_dict  in model_params_prev:
                unique_param = model_param_dict['param_name']+'.'+\
                                    model_param_dict['sectionlist']
                model_param_opt_entry = list(filter(lambda x: x['param_name']+'.'+\
                                    x['sectionlist'] == unique_param, model_params_opt))[0]

                model_param_opt_entry = adjust_param_rule(model_param_opt_entry,model_param_dict)


        # Add reversal potential if Na, K currents are present

        rev_potential =  utility.rev_potential
        for rev in rev_potential:
            if rev == 'ena':
                for sect in ena_sect:
                    if self.no_apical and sect == 'apic': # if no apical dendrite in morphology
                        continue
                    iter_dict =  {'param_name':rev, 'sectionlist':section_map[sect], 'dist_type': 'uniform',
                              'type':'section','value':rev_potential[rev]}
                    model_params_opt.append(iter_dict)
            elif rev == 'ek':
                for sect in ek_sect:
                    if self.no_apical and sect == 'apic': # if no apical dendrite in morphology
                        continue
                    iter_dict = {'param_name':rev, 'sectionlist':section_map[sect], 'dist_type': 'uniform',
                              'type':'section','value':rev_potential[rev]}
                    model_params_opt.append(iter_dict)

        # Add experimental conditions

        model_params_opt.append({"param_name": "celsius","type": "global","value": self.temperature})
        model_params_opt.append({"param_name": "v_init","type": "global","value": self.v_init})
        model_params_release = self.get_release_params(section_map,rev_potential)

        return model_params_opt, model_params_release



    def load_params_prev_stage(self,section_map):
        model_prev_stage = utility.load_json(self.prev_stage_model_path)
        model_params_prev = list()
        for key, values in model_prev_stage.items():
            if key == 'genome':
                for j in range(len(values)):
                    if self.no_apical and model_prev_stage[key][j]['section'] == 'apic': # if no apical dendrite in morphology
                        continue

                    iter_dict = {'param_name':model_prev_stage[key][j]['name']}
                    iter_dict['dist_type'] = 'uniform'
                    iter_dict['sectionlist'] = section_map[model_prev_stage[key][j]['section']]
                    iter_dict['value'] = float(model_prev_stage[key][j]['value'])
                    iter_dict['type'] = 'section'
                    if model_prev_stage[key][j]['mechanism'] != '':
                        iter_dict['mech'] = model_prev_stage[key][j]['mechanism']
                        iter_dict['type'] = 'range'
                    model_params_prev.append(iter_dict)
        return model_params_prev


    def get_release_params(self,section_map,rev_potential):

        if self.release_param_path:
            model_params_release = list()
            data_release = utility.load_json(self.release_param_path)
            for key, values in data_release.items():
                if key == 'genome':
                    for j in range(len(values)):
                             iter_dict_release = {'param_name':data_release[key][j]['name']}
                             iter_dict_release['sectionlist'] = section_map[data_release[key][j]['section']]
                             iter_dict_release['type'] = 'section'
                             iter_dict_release['value'] = float(data_release[key][j]['value'])
                             iter_dict_release['dist_type'] = 'uniform'
                             if data_release[key][j]['mechanism'] != '':
                                    iter_dict_release['mech'] = data_release[key][j]['mechanism']
                                    iter_dict_release['type'] = 'range'
                             model_params_release.append(iter_dict_release)

            for sect in list(set(section_map.values())-set(['all'])):
                for rev in rev_potential:
                    iter_dict_release =  {'param_name':rev, 'sectionlist':sect, 'dist_type': 'uniform', 'type':'section'}
                    if rev == 'ena':
                        iter_dict_release['value'] = rev_potential[rev]
                    elif rev == 'ek':
                        iter_dict_release['value'] = rev_potential[rev]
                    model_params_release.append(iter_dict_release)
            model_params_release.append({"param_name": "celsius","type": "global","value": self.temperature})
            model_params_release.append({"param_name": "v_init","type": "global","value": -90})
        else:
            model_params_release = None

        return model_params_release

    def write_params_opt(self,model_params,model_params_release,
                         base_dir = 'config/',**kwargs):

        param_write_path =  kwargs.get('param_write_path') or \
                base_dir+ self.cell_id + '/parameters.json'
        release_param_write_path = kwargs.get('release_param_write_path') or \
            base_dir+ self.cell_id + '/release_parameters.json'

        utility.create_filepath(param_write_path)
        utility.save_json(param_write_path,model_params)

        release_params = dict() # for comparison with optimized values (strictly for values)
        if model_params_release:
            for param_dict_release in model_params_release:
                param_name = param_dict_release['param_name']
                if param_name not in ['ena','ek','v_init','celsius']:
                    release_params[param_name + '.' + param_dict_release['sectionlist']] = param_dict_release['value']

            # Released parameter file in bpopt format for running simulation
            utility.save_json(release_param_write_path,model_params_release)
        else:
            release_param_write_path = None

        return param_write_path,release_param_write_path,release_params

    def get_opt_mechanism(self,model_params,model_params_release,
                             param_bounds_path):
        params_dict = utility.load_json(param_bounds_path)
        active_params, Ih_params, _,_,_,_=self.group_params(params_dict)
        model_mechs = defaultdict(list)
        model_mechs['all'].append('pas')

        for param_dict in model_params:
            if param_dict['param_name'] in active_params+Ih_params:
                if param_dict['mech'] not in model_mechs[param_dict['sectionlist']]:
                    model_mechs[param_dict['sectionlist']].append(param_dict['mech'])

        if model_params_release:
            model_mechs_release = {'somatic' : ['pas'], 'axonal':['pas'], 'apical':['pas'],
                               'basal': ['pas']}
            for param_dict_release in model_params_release:
                if 'mech' in param_dict_release.keys():
                    if param_dict_release['mech'] not in model_mechs_release[param_dict_release['sectionlist']]:
                        model_mechs_release[param_dict_release['sectionlist']].append(param_dict_release['mech'])
        else:
            model_mechs_release = None

        return model_mechs,model_mechs_release

    def write_mechanisms_opt(self,model_mechs,model_mechs_release,
                             base_dir = 'config/',**kwargs):

        mechanism_write_path =  kwargs.get('mechanism_write_path') or \
            base_dir+ self.cell_id + '/mechanism.json'
        mechanism_release_write_path = kwargs.get('mechanism_release_write_path') \
            or base_dir+ self.cell_id + '/mechanism_release.json'

        utility.create_filepath(mechanism_write_path)
        utility.save_json(mechanism_write_path,model_mechs)

        if model_mechs_release:
            utility.save_json(mechanism_release_write_path,model_mechs_release)
        else:
            model_mechs_release = None
            mechanism_release_write_path = None

        return mechanism_write_path,mechanism_release_write_path

    def aibs_peri_to_bpopt(self,peri_param_path,base_dir='config/'):
        peri_params = utility.load_json(peri_param_path)
        peri_params_release = list()
        peri_mechs_release = defaultdict(list)
        peri_mechs_release['all'].append('pas')

        rev_potential =  utility.rev_potential
        section_map = utility.bpopt_section_map

        for key, values in peri_params.items():
            if key == 'genome':
                for j in range(len(values)):
                     iter_dict_release = {'param_name':peri_params[key][j]['name']}
                     iter_dict_release['sectionlist'] = section_map[peri_params[key][j]['section']]
                     iter_dict_release['type'] = 'section'
                     iter_dict_release['value'] = float(peri_params[key][j]['value'])
                     iter_dict_release['dist_type'] = 'uniform'
                     if peri_params[key][j]['mechanism'] != '':
                            iter_dict_release['mech'] = peri_params[key][j]['mechanism']
                            iter_dict_release['type'] = 'range'
                     peri_params_release.append(iter_dict_release)

            elif key == 'passive':
                 for key_pas,val_pas in values[0].items():
                     if key_pas == 'cm':
                         for pas_param in val_pas:
                             iter_dict_release ={'param_name':'cm',
                                                 'sectionlist' : section_map[pas_param['section']],
                                                 'value' : pas_param['cm'],
                                                 'dist_type' : 'uniform',
                                                 'type' : 'section'}
                             peri_params_release.append(iter_dict_release)

                     else:
                          iter_dict_release ={'param_name':'Ra' \
                                              if key_pas== 'ra' else key_pas,
                                                 'sectionlist' : 'all',
                                                 'value' : val_pas,
                                                 'dist_type' : 'uniform',
                                                 'type' : 'section'}
                          peri_params_release.append(iter_dict_release)


        for rev in rev_potential:
            iter_dict_release =  {'param_name':rev, 'sectionlist':'somatic',
                                  'dist_type': 'uniform', 'type':'section'}
            if rev == 'ena':
                iter_dict_release['value'] = rev_potential[rev]
            elif rev == 'ek':
                iter_dict_release['value'] = rev_potential[rev]
            peri_params_release.append(iter_dict_release)

        peri_params_release.append({"param_name": "celsius","type": "global","value": 34})
        peri_params_release.append({"param_name": "v_init","type": "global",
                                     "value": peri_params['conditions'][0]["v_init"]})

        for param_dict in peri_params_release:
            if 'mech' in param_dict.keys():
                if param_dict['mech'] not in peri_mechs_release[param_dict['sectionlist']]:
                    peri_mechs_release[param_dict['sectionlist']].append(param_dict['mech'])

        peri_params_write_path = base_dir+ self.cell_id + '/peri_parameters.json'
        peri_mech_write_path = base_dir+ self.cell_id + '/peri_mechanism.json'
        utility.save_json(peri_params_write_path,peri_params_release)
        utility.save_json(peri_mech_write_path,peri_mechs_release)
        return peri_params_write_path, peri_mech_write_path


    def write_opt_config_file(self,morph_path,param_write_path,
                              mech_write_path,mech_release_write_path,
                              features_write_path,untrained_features_write_path,
                              all_features_write_path,
                              protocols_write_path,all_protocols_write_path,
                              release_params,release_param_write_path,
                              opt_config_filename = 'config_file.json',
                              **kwargs):
        """Writes a bunch of paths to json. Can also hack to write other keys for the Evaluator"""

        path_dict =  dict()
        path_dict['morphology'] = morph_path
        path_dict['parameters'] = param_write_path
        path_dict['mechanism'] = mech_write_path
        path_dict['mechanism_release'] = mech_release_write_path
        path_dict['features'] = features_write_path
        path_dict['untrained_features'] = untrained_features_write_path
        path_dict['all_features'] = all_features_write_path
        path_dict['protocols'] = protocols_write_path
        path_dict['all_protocols'] = all_protocols_write_path
        path_dict['release_params_bpopt'] = release_params
        path_dict['released_model'] = release_param_write_path
        for config_key,path in kwargs.items():
            path_dict[config_key] = path

        utility.save_json(opt_config_filename,path_dict)
