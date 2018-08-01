#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 13:48:46 2018

@author: anin
"""

import json

import model_helper  # NOQA

import bluepyopt.ephys as ephys
import logging
logger = logging.getLogger(__name__)



# TODO store definition dicts in json
# TODO rename 'score' into 'objective'
# TODO add functionality to read settings of every object from config format


def define_protocols(protocols_write_path):
    """Define protocols"""

    protocol_definitions = json.load(open(protocols_write_path))

    protocols = {}

    soma_loc = ephys.locations.NrnSeclistCompLocation(
        name='soma',
        seclist_name='somatic',
        sec_index=0,
        comp_x=0.5)
    
    AIS_loc = ephys.locations.NrnSeclistCompLocation(
        name='AIS',
        seclist_name='axonal',
        sec_index=1,
        comp_x=0.5)
    
    for protocol_name, protocol_definition in protocol_definitions.items():
        
        # By default include somatic and AIS recording
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
        
        recordings = [somav_recording, AIS_recording]

        if 'extra_recordings' in protocol_definition:
            for recording_definition in protocol_definition['extra_recordings']:
                if recording_definition['type'] == 'somadistance':
                    location = ephys.locations.NrnSomaDistanceCompLocation(
                        name=recording_definition['name'],
                        soma_distance=recording_definition['somadistance'],
                        seclist_name=recording_definition['seclist_name'])
                    var = recording_definition['var']
                    recording = ephys.recordings.CompRecording(
                        name='%s.%s.%s' % (protocol_name, location.name, var),
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
                    ramp_amplitude_end = stimulus_definition['amp_end'],
                    ramp_delay=stimulus_definition['delay'],
                    ramp_duration=stimulus_definition['duration'],
                    location=soma_loc,
                    total_duration=stimulus_definition['totduration']))
                    

        protocols[protocol_name] = ephys.protocols.SweepProtocol(
            protocol_name,
            stimuli,
            recordings)

    return protocols


def define_fitness_calculator(protocols, features_write_path, weight_dict):
    """Define fitness calculator"""

    feature_definitions = json.load(open(features_write_path))

    # TODO: add bAP stimulus
    objectives = []

    for protocol_name, locations in feature_definitions.items():
        for location, features in locations.items():
#            features_vec = []
#            weighted_feature_name = protocol_name + '.' + location
#            weight_vec = list()
            
            for efel_feature_name, meanstd in features.items():
                feature_name = '%s.%s.%s' % (protocol_name, location, efel_feature_name)
#                weighted_feature_name += '.' + efel_feature_name 
                recording_names = {'': '%s.%s.v' % (protocol_name, location), 
                                   'location_AIS' : '%s.AIS.v' %protocol_name}
                stimulus = protocols[protocol_name].stimuli[0]
                if 'Ramp' in protocol_name:
                    stim_start = stimulus.ramp_delay 
                    duration = stimulus.ramp_duration
                elif 'DC' in protocol_name:
                    stim_start = stimulus.step_delay
                    duration = stimulus.step_duration

                if location == 'soma':
                    threshold = -20
                elif 'dend' in location:
                    threshold = -55

                if protocol_name == 'bAP':
                    stim_end = stimulus.total_duration
#                    stim_end = stimulus.step_delay + stimulus.step_duration

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
#                if weight_dict:
#                    weight_vec.append(weight_dict[efel_feature_name])
#                    features_vec.append(feature)
				
#                else:
                objective = ephys.objectives.SingletonObjective(
                    feature_name,
                    feature)
                objectives.append(objective)
#            if weight_dict:
#                objective = ephys.objectives.WeightedSumObjective(weighted_feature_name,
#                                                                  features_vec,weight_vec)
#                objectives.append(objective)

    fitcalc = ephys.objectivescalculators.ObjectivesCalculator(objectives)

    return fitcalc


def create(protocol_path,feature_path, morph_path, param_path, mech_path, weight_dict = None):
    """Setup"""

    cell = model_helper.create(morph_path, param_path,mech_path)

    fitness_protocols = define_protocols(protocol_path)
    fitness_calculator = define_fitness_calculator(fitness_protocols,feature_path, weight_dict)

    param_names = [param.name
                   for param in cell.params.values()
                   if not param.frozen]
    
    sim = ephys.simulators.NrnSimulator()

    return ephys.evaluators.CellEvaluator(
        cell_model=cell,
        param_names=param_names,
        fitness_protocols=fitness_protocols,
        fitness_calculator=fitness_calculator,
        sim=sim)
