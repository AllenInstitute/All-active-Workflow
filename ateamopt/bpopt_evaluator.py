import json
import numpy as np
import bluepyopt.ephys as ephys
from ateamopt.utils import utility
import logging
import os

logger = logging.getLogger(__name__)


class Bpopt_Evaluator(object):

    def __init__(self, protocol_path, feature_path,
                 morph_path, param_path, mech_path, ephys_dir='preprocessed',
                 skip_features=['peak_time'],
                 **props):
        """
        do_replace_axon : bluepyopt axon replace code, diameter taken
        from swc file
        """
        self.morph_path = morph_path if morph_path else None
        self.protocol_path = protocol_path if protocol_path else None
        self.feature_path = feature_path if feature_path else None
        self.param_path = param_path if param_path else None
        self.mech_path = mech_path if mech_path else None
        self.ephys_dir = ephys_dir
        self.AIS_check = False
        self.skip_features = skip_features

        if self.feature_path:
            feature_definitions = utility.load_json(self.feature_path)
            feature_set = []
            for feat_key, feat_val in feature_definitions.items():
                feature_set.extend(feat_val['soma'].keys())
            self.AIS_check = True if 'check_AISInitiation' in \
                list(set(feature_set)) else False
        if any(timed_prop in props for timed_prop in ['timeout', 'learn_eval_trend']):
            self.timed_evaluation = True
        else:
            self.timed_evaluation = False

        self.axon_type = 'stub_axon'
        if props.pop('do_replace_axon', None):
            self.axon_type = 'bpopt_replaced_axon'

        self.eval_props = props

    def define_mechanisms(self):
        """Define mechanisms"""

        mech_path = self.mech_path

        mech_definitions = json.load(
            open(mech_path))

        mechanisms = []
        for sectionlist, channels in mech_definitions.items():
            seclist_loc = ephys.locations.NrnSeclistLocation(
                sectionlist,
                seclist_name=sectionlist)
            for channel in channels:
                mechanisms.append(ephys.mechanisms.NrnMODMechanism(
                    name='%s.%s' % (channel, sectionlist),
                    mod_path=None,
                    suffix=channel,
                    locations=[seclist_loc],
                    preloaded=True))

        return mechanisms

    def define_parameters(self):
        """Define parameters"""

        param_path = self.param_path
        param_configs = json.load(open(param_path))
        parameters = []

        for param_config in param_configs:
            if 'value' in param_config:
                frozen = True
                value = param_config['value']
                bounds = None
            elif 'bounds' in param_config:
                frozen = False
                bounds = param_config['bounds']
                value = None
            else:
                raise Exception(
                    'Parameter config has to have bounds or value: %s'
                    % param_config)

            if param_config['type'] == 'global':
                parameters.append(
                    ephys.parameters.NrnGlobalParameter(
                        name=param_config['param_name'],
                        param_name=param_config['param_name'],
                        frozen=frozen,
                        bounds=bounds,
                        value=value))
            elif param_config['type'] in ['section', 'range']:
                if param_config['dist_type'] == 'uniform':
                    scaler = ephys.parameterscalers.NrnSegmentLinearScaler()
                elif param_config['dist_type'] == 'exp':
                    scaler = ephys.parameterscalers.NrnSegmentSomaDistanceScaler(
                        distribution=param_config['dist'])
                seclist_loc = ephys.locations.NrnSeclistLocation(
                    param_config['sectionlist'],
                    seclist_name=param_config['sectionlist'])

                name = '%s.%s' % (param_config['param_name'],
                                  param_config['sectionlist'])

                if param_config['type'] == 'section':
                    parameters.append(
                        ephys.parameters.NrnSectionParameter(
                            name=name,
                            param_name=param_config['param_name'],
                            value_scaler=scaler,
                            value=value,
                            frozen=frozen,
                            bounds=bounds,
                            locations=[seclist_loc]))
                elif param_config['type'] == 'range':
                    parameters.append(
                        ephys.parameters.NrnRangeParameter(
                            name=name,
                            param_name=param_config['param_name'],
                            value_scaler=scaler,
                            value=value,
                            frozen=frozen,
                            bounds=bounds,
                            locations=[seclist_loc]))
            else:
                raise Exception(
                    'Param config type has to be global, section or range: %s' %
                    param_config)

        return parameters

    def define_morphology(self):
        """Define morphology"""

        morph_path = self.morph_path
        return ephys.morphologies.NrnFileMorphology(morph_path,
                                                    stub_axon=self.axon_type == 'stub_axon',
                                                    do_replace_axon=self.axon_type != 'stub_axon')

    def model_builder(self):
        """Create cell model"""

        cell = ephys.models.CellModel(
            'cell',
            morph=self.define_morphology(),
            mechs=self.define_mechanisms(),
            params=self.define_parameters())

        return cell

    def define_protocols(self):
        """Define protocols, loading definitions from self.protocol_path
        Returns: 
            protocols (dict of str: bluepyopt.ephys.protocols.SweepProtocol):
                protocols indexed by name
        """
        ephys_dir = self.ephys_dir
        protocol_definitions = json.load(open(self.protocol_path))

        protocols = {}

        soma_loc = ephys.locations.NrnSeclistCompLocation(
            name='soma',
            seclist_name='somatic',
            sec_index=0,
            comp_x=0.5)

        AIS_loc = ephys.locations.NrnSeclistCompLocation(
            name='AIS',
            seclist_name='axonal',
            sec_index=0,
            comp_x=0.5)

        for protocol_name, protocol_definition in protocol_definitions.items():
            # By default include somatic recording
            somav_recording = ephys.recordings.CompRecording(
                name='%s.soma.v' %
                protocol_name,
                location=soma_loc,
                variable='v')

            AIS_recording = ephys.recordings.CompRecording(
                name='%s.AIS.v' %
                protocol_name,
                location=AIS_loc,
                variable='v')

            if self.AIS_check and 'LongDC' in protocol_name: 
                # Only for Square pulses (AIS recording removed for others later)
                # TODO: Possibly pass this through config file(which stimtypes should record at AIS)
                recordings = [somav_recording, AIS_recording]
            else:
                recordings = [somav_recording]

            if 'extra_recordings' in protocol_definition:
                for recording_definition in protocol_definition['extra_recordings']:
                    if recording_definition['type'] == 'somadistance':
                        location = ephys.locations.NrnSomaDistanceCompLocation(
                            name=recording_definition['name'],
                            soma_distance=recording_definition['somadistance'],
                            seclist_name=recording_definition['seclist_name'])
                        var = recording_definition['var']
                        recording = ephys.recordings.CompRecording(
                            name='%s.%s.%s' % (
                                protocol_name, location.name, var),
                            location=location,
                            variable=recording_definition['var'])

                        recordings.append(recording)
                    else:
                        raise Exception(
                            'Recording type %s not supported' %
                            recording_definition['type'])

            stimuli = []
            for stimulus_definition in protocol_definition['stimuli']:
                if stimulus_definition['type'] == 'SquarePulse':
                    stimuli.append(ephys.stimuli.NrnSquarePulse(
                        step_amplitude=stimulus_definition['amp'],
                        step_delay=stimulus_definition['delay'],
                        step_duration=stimulus_definition['duration'],
                        location=soma_loc,
                        total_duration=stimulus_definition['totduration']))
                elif stimulus_definition['type'] == 'RampPulse':
                    stimuli.append(ephys.stimuli.NrnRampPulse(
                        ramp_amplitude_start=stimulus_definition['amp'],
                        ramp_amplitude_end=stimulus_definition['amp_end'],
                        ramp_delay=stimulus_definition['delay'],
                        ramp_duration=stimulus_definition['duration'],
                        location=soma_loc,
                        total_duration=stimulus_definition['totduration']))
                    
                elif stimulus_definition['type'] in ['TriBlip', 'Noise']:
                    sweep_file = os.path.join(ephys_dir,
                                              stimulus_definition['sweep_filenames'][0])
                    stim_play_time = np.loadtxt(sweep_file)[:, 0]
                    stim_play_current = np.loadtxt(sweep_file)[:, 2]
                    stimuli.append(ephys.stimuli.NrnCurrentPlayStimulus(
                        current_points=stim_play_current,
                        time_points=stim_play_time,
                        location=soma_loc))
                    
                   
            protocols[protocol_name] = ephys.protocols.SweepProtocol(
                protocol_name,
                stimuli,
                recordings)

        return protocols

    def define_fitness_calculator(self, fitness_protocols):
        """Define fitness calculator, loading definitions from self.feature_path"""

        # TODO: add bAP stimulus
        objectives = []
        feature_definitions = json.load(open(self.feature_path))
        protocol_definitions = json.load(open(self.protocol_path))
        for protocol_name, locations in feature_definitions.items():
            for location, features in locations.items():

                for efel_feature_name, meanstd in features.items():
                    if efel_feature_name not in self.skip_features:
                        feature_name = '%s.%s.%s' % (
                            protocol_name, location, efel_feature_name)
                        if self.AIS_check and 'LongDC' in protocol_name:
                            recording_names = {'': '%s.%s.v' % (protocol_name, location),
                                               'location_AIS': '%s.AIS.v' % protocol_name}
                        else:
                            recording_names = {
                                '': '%s.%s.v' % (protocol_name, location)}
                        stimulus = fitness_protocols[protocol_name].stimuli[0]
                        if 'Ramp' in protocol_name:
                            stim_start = stimulus.ramp_delay
                            duration = stimulus.ramp_duration
                        elif 'DC' in protocol_name:
                            stim_start = stimulus.step_delay
                            duration = stimulus.step_duration
                        else:
                            stim_start = protocol_definitions[protocol_name]['stimuli'][0]['delay']
                            duration = protocol_definitions[protocol_name]['stimuli'][0]['duration']
                        
                        if location == 'soma':
                            threshold = -20
                        elif 'dend' in location:
                            threshold = -55
    
                        if protocol_name == 'bAP':
#                            stim_end = stimulus.total_duration
                            stim_end = stimulus.step_delay + stimulus.step_duration
    
                        else:
                            stim_end = stim_start + duration
    
                        feature = ephys.efeatures.eFELFeature(
                            feature_name,
                            efel_feature_name=efel_feature_name,
                            recording_names=recording_names,
                            stim_start=stim_start,
                            stim_end=stim_end,
                            exp_mean=meanstd[0],
                            exp_std=meanstd[1],
                            threshold=threshold,
                            force_max_score=True,
                            max_score=250)
    
                        objective = ephys.objectives.SingletonObjective(
                            feature_name,
                            feature)
                        objectives.append(objective)

        fitcalc = ephys.objectivescalculators.ObjectivesCalculator(objectives)

        return fitcalc

    def create_evaluator(self):
        """Setup"""

        cell = self.model_builder()

        fitness_protocols = self.define_protocols()
        fitness_calculator = self.define_fitness_calculator(fitness_protocols)

        param_names = [param.name
                       for param in cell.params.values()
                       if not param.frozen]

        sim = ephys.simulators.NrnSimulator()

        if self.timed_evaluation:
            kwargs = {}
            for key, val in self.eval_props.items():
                if val:
                    kwargs[key] = val

            return ephys.evaluators.CellEvaluatorTimed(
                cell_model=cell,
                param_names=param_names,
                fitness_protocols=fitness_protocols,
                fitness_calculator=fitness_calculator,
                sim=sim, **kwargs)
        else:
            return ephys.evaluators.CellEvaluator(
                cell_model=cell,
                param_names=param_names,
                fitness_protocols=fitness_protocols,
                fitness_calculator=fitness_calculator,
                sim=sim)
