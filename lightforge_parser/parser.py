#
# This is a parser for Nanomatch GmbH's kMC software "Lightforge"
#
#
#
#
#
import yaml
import os
import re
import datetime
import numpy as np
from pathlib import Path

from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg as units
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.system import System
from nomad.datamodel.metainfo.simulation.calculation import Calculation, Energy, EnergyEntry
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.datamodel.results import Results, Properties, Structure
from nomad.parsing.file_parser import UnstructuredTextFileParser, Quantity
from nomad.datamodel.optimade import Species
from . import metainfo  # pylint: disable=unused-import
from .metainfo.lightforge import (
                            IV, IQE2, Current_density, Current_characteristics, Experiments, Material,
                            Input, Mobility, Particle_densities, Charge_density_average, 
                            Exciton_decay_density_average, Photon_creation_density_average,
                            Quenching_density_average, Exciton_molpairs, Emitter_emitter_transport_count,
                            Host_emitter_transport_count, Host_host_transport_count, Runtime_analysis, Event_counts_by_type, Device_data, Electrodes, Energy_levels,
                            Exciton_separation, Foerster, Site_energies, Mol_types, Coordinates) 


def DetailedParser(filepath, archive):
    sec_run = archive.m_create(Run)
    sec_calc = sec_run.m_create(Calculation)
    sec_experiments =  sec_calc.m_create(Experiments)
    sec_material = sec_calc.m_create(Material)
    sec_current_characteristics = sec_experiments.m_create(Current_characteristics)
    sec_particle_densities = sec_experiments.m_create(Particle_densities)
    
    sec_IQE2 = sec_current_characteristics.m_create(IQE2)
    sec_IV = sec_current_characteristics.m_create(IV)
    exciton_molpairs_hasrun = False
    runtime_analysis_hasrun = False
    device_data_hasrun= False
    foerster_hasrun = False
    for root, dirs, files in sorted(os.walk(filepath.parent)):
        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]
        
        files = sorted(files, key = natsort)
        i = 0
        while i < len(files):
            if '.png' in files[i] or '.npz' in files[i]:
                files.remove(files[i])
            else:
                i += 1    
        for file in files:    
            with open(root +'/'+ file, 'r') as f:
                if 'current_density' in file and 'all_data_points' not in root:
                    sec_current_density = sec_current_characteristics.m_create(Current_density)
                    value = []
                    for i, line in enumerate(f):
                        line = float(line)
                        value.append(line)
                    sec_current_density.value = np.array(value)
                if 'IQE2_all_currents' in file and 'all_data_points' not in root:
                    for i, line in enumerate(f):
                        rows = i+1
                    
                    a = np.zeros((rows,2))
                    sec_IQE2.iqe2_all_currents = a

                    f.seek(0)
                    for i, line in enumerate(f):
                        
                        parts = line.split()
                        
                        sec_IQE2.iqe2_all_currents[i][0] = parts[0]
                        sec_IQE2.iqe2_all_currents[i][1] = parts[1]
                if 'IQE2_all_fields' in file and 'all_data_points' not in root:
                    for i, line in enumerate(f):
                        rows = i + 1
                    b = np.zeros((rows,2))
                    sec_IQE2.iqe2_all_fields = b
                    f.seek(0)
                    for i, line in enumerate(f):
                        parts = line.split()
                        sec_IQE2.iqe2_all_fields[i][0] = parts[0]
                        sec_IQE2.iqe2_all_fields[i][1] = parts[1] 
                if re.search(r'^IV_all_fields.dat$', file) and 'all_data_points' not in root:
                    for i, line in enumerate(f):
                        rows = i + 1
                    c = np.zeros((rows,3))
                    sec_IV.iv_all_fields = c
                    f.seek(0)
                    for i, line in enumerate(f):
                        parts = line.split()
                        
                        sec_IV.iv_all_fields[i][0] = parts[0]
                        sec_IV.iv_all_fields[i][1] = parts[1]                 
                        sec_IV.iv_all_fields[i][2] = parts[2]           
                if re.search(r'mobilities_\d+', file) and 'all_data_points' not in root:
                    sec_mobility = sec_current_characteristics.m_create(Mobility)
                    value = []
                    for i, line in enumerate(f):
                        line=float(line)
                        value.append(line)
                    sec_mobility.value = np.array(value)
                if re.search(r'mobilities_all_fields', file) and 'all_data_points' not in root:
                    sec_mobility = sec_current_characteristics.m_create(Mobility)
                    for i, line in enumerate(f):
                        rows  = i +1 
                    d = np.zeros((rows, 3))
                    sec_mobility.mobilities_all_fields = d
                    f.seek(0)
                    for i, line in enumerate(f):
                        parts = line.split()
                        sec_mobility.mobilities_all_fields[i][0] = parts[0]
                        sec_mobility.mobilities_all_fields[i][1] = parts[1]
                        sec_mobility.mobilities_all_fields[i][2] = parts[2]
                if re.search(r'charge_density_average_\d+.dat', file) and 'all_data_points' not in root:
                    sec_charge_density_average = sec_particle_densities.m_create(Charge_density_average)
                    device_length = []
                    electrons = []
                    holes = []
                    for i, line in enumerate(f):
                        columns = len(line.split())
                        break
                    f.seek(0)
                    value = np.zeros((3, columns))
                    for i, line in enumerate(f):
                        if i==0:
                            parts = line.split()
                            device_length = parts
                            
                        if i==1:
                            parts = line.split()
                            electrons = parts
                            
                        if i==2:
                            parts = line.split()
                            holes = parts
                    value[0] = device_length
                    value[1] = electrons
                    value[2] = holes  
                    sec_charge_density_average.value = value
                
                if re.search(r'exciton_decay_density_average_\d+', file)  and 'all_data_points'  not in root:
                    sec_exciton_decay_density_average  = sec_particle_densities.m_create(Exciton_decay_density_average)
                    file_exciton_decay_density_average =  yaml.safe_load(f)
                                        
                    if 'CS' in file_exciton_decay_density_average['total']['annihilation']:
                        _cs = file_exciton_decay_density_average['total']['annihilation']['CS']
                        sec_exciton_decay_density_average.cs = _cs    
                    if 'EET' in file_exciton_decay_density_average['total']['annihilation']:
                        _eet = file_exciton_decay_density_average['total']['annihilation']['EET']
                        sec_exciton_decay_density_average.eet = _eet
                    if 'EPT' in file_exciton_decay_density_average['total']['annihilation']:
                        _ept = file_exciton_decay_density_average['total']['annihilation']['EPT']
                        sec_exciton_decay_density_average.ept = _ept                    
                    if 'PTQ' in file_exciton_decay_density_average['total']['annihilation']:
                        _ptq = file_exciton_decay_density_average['total']['annihilation']['PTQ']
                        sec_exciton_decay_density_average.ptq = _ptq
                    if 'SPQ' in file_exciton_decay_density_average['total']['annihilation']:
                        _spq = file_exciton_decay_density_average['total']['annihilation']['SPQ']
                        sec_exciton_decay_density_average.spq= _spq
                    if 'SSA' in file_exciton_decay_density_average['total']['annihilation']:
                        _ssa = file_exciton_decay_density_average['total']['annihilation']['SSA']
                        sec_exciton_decay_density_average.ssa = _ssa
                    if 'STA' in file_exciton_decay_density_average['total']['annihilation']:
                        _sta = file_exciton_decay_density_average['total']['annihilation']['STA']
                        sec_exciton_decay_density_average.sta = _sta
                    if 'TPQ' in file_exciton_decay_density_average['total']['annihilation']:
                        _tpq = file_exciton_decay_density_average['total']['annihilation']['TPQ']
                        sec_exciton_decay_density_average.tpq = _tpq
                    if 'TSA' in file_exciton_decay_density_average['total']['annihilation']:
                        _tsa = file_exciton_decay_density_average['total']['annihilation']['TSA']
                        sec_exciton_decay_density_average.tsa = _tsa 
                    if 'TTA' in file_exciton_decay_density_average['total']['annihilation']:
                        _tta = file_exciton_decay_density_average['total']['annihilation']['TTA']
                        sec_exciton_decay_density_average.tta = _tta
                    if 'TTF' in file_exciton_decay_density_average['total']['annihilation']:
                        _ttf = file_exciton_decay_density_average['total']['annihilation']['TTF']
                        sec_exciton_decay_density_average.ttf = _ttf
                    if 'radiative' in file_exciton_decay_density_average['total']['annihilation']:
                        _radiative = file_exciton_decay_density_average['total']['annihilation']['radiative']
                        sec_exciton_decay_density_average.radiative = _radiative
                    if 'thermal' in file_exciton_decay_density_average['total']['annihilation']:
                        _thermal = file_exciton_decay_density_average['total']['annihilation']['thermal']
                        sec_exciton_decay_density_average.thermal = _thermal
                    if 'photon' in file_exciton_decay_density_average['total']['creation']:
                        _photon = file_exciton_decay_density_average['total']['creation']['photon']
                        sec_exciton_decay_density_average.photon = _photon
                    if 'recombination' in file_exciton_decay_density_average['total']['creation']:
                        _recombination = file_exciton_decay_density_average['total']['creation']['recombination']
                        sec_exciton_decay_density_average.recombination = _recombination
                    if 'x_axis' in file_exciton_decay_density_average:
                        _x_axis = file_exciton_decay_density_average['x_axis']
                        sec_exciton_decay_density_average.x_axis = _x_axis    
                
                if re.search(r'photon_creation_density_average_\d+', file) and 'all_data_points' not in root:
                    sec_photon_creation_density_average = sec_particle_densities.m_create(Photon_creation_density_average)
                    device_length = []
                    photons = []
                    excitons = []
                    for i, line in enumerate(f):
                        columns = len(line.split())
                        break
                    f.seek(0)
                    value = np.zeros((3, columns))
                    for i, line in enumerate(f):
                        if i == 0:
                            parts = line.split()
                            device_length = parts
                        if i == 1:
                            parts = line.split()
                            photons = parts
                        if i == 2:
                            parts = line.split()
                            excitons = parts
                    value[0] = device_length
                    value[1] = photons
                    value[2] = excitons
                    sec_photon_creation_density_average.value = value
                if re.search(r'quenching_density_average_\d+', file) and 'all_data_points' not in root:
                    sec_quenching_density_average = sec_particle_densities.m_create(Quenching_density_average)
                    device_length = []
                    excitons_quenched = []
                    electrons = []
                    holes = []
                    for i, line in enumerate(f):
                        columns = len(line.split())
                        break
                    f.seek(0)
                    value = np.zeros((4, columns))
                    for i, line in enumerate(f):
                        if i == 0:
                            parts = line.split()
                            device_length = parts
                        if i == 1:
                            parts = line.split()
                            excitons_quenched = parts
                        if i == 2:
                            parts = line.split()
                            electrons = parts
                        if i == 3:
                            parts = line.split()
                            holes = parts
                    value[0] = device_length     
                    value[1] = excitons_quenched
                    value[2] = electrons
                    value[3] = holes 
                    sec_quenching_density_average.value = value
                if 'exciton_molpairs' in root:
                    if not exciton_molpairs_hasrun:
                        sec_exciton_molpairs = sec_particle_densities.m_create(Exciton_molpairs)
                        exciton_molpairs_hasrun = True   
                    if re.search(r'host_emitter_transport_count', file) and 'all_data_points' not in root:
                        sec_host_emitter_transport_count = sec_exciton_molpairs.m_create(Host_emitter_transport_count) # might be wrong
                        file_host_emitter_transport_count = yaml.safe_load(f)
                        if 'dexter_S1S1' in file_host_emitter_transport_count:
                            _dexter_s1s1 = file_host_emitter_transport_count['dexter_S1S1']
                            sec_host_emitter_transport_count.dexter_s1s1 = _dexter_s1s1
                        if 'dexter_S1T1' in file_host_emitter_transport_count:
                            _dexter_s1t1 = file_host_emitter_transport_count['dexter_S1T1']
                            sec_host_emitter_transport_count.dexter_s1t1 = _dexter_s1t1
                        if 'dexter_T1S1' in file_host_emitter_transport_count:
                            _dexter_t1s1 = file_host_emitter_transport_count['dexter_T1S1']
                            sec_host_emitter_transport_count.dexter_t1s1 = _dexter_t1s1
                        if 'dexter_T1T1' in file_host_emitter_transport_count:
                            _dexter_t1t1 = file_host_emitter_transport_count['dexter_T1T1']
                            sec_host_emitter_transport_count.dexter_t1t1 = _dexter_t1t1
                        if 'foerster_S1S1' in file_host_emitter_transport_count:
                            _foerster_s1s1 = file_host_emitter_transport_count['foerster_S1S1']
                            sec_host_emitter_transport_count.foerster_s1s1 = _foerster_s1s1
                        if 'foerster_S1T1' in file_host_emitter_transport_count:
                            _foerster_s1t1 = file_host_emitter_transport_count['foerster_S1T1']
                            sec_host_emitter_transport_count.foerster_s1t1 = _foerster_s1t1
                        if 'foerster_T1S1' in file_host_emitter_transport_count:
                            _foerster_t1s1 = file_host_emitter_transport_count['foerster_T1S1']
                            sec_host_emitter_transport_count.foerster_t1s1 = _foerster_t1s1
                        if 'foerster_T1T1' in file_host_emitter_transport_count:
                            _foerster_t1t1 = file_host_emitter_transport_count['foerster_T1T1']
                            sec_host_emitter_transport_count.foerster_t1t1 = _foerster_t1t1
                        if 'x_axis' in file_host_emitter_transport_count:
                            _x_axis = file_host_emitter_transport_count['x_axis']
                            sec_host_emitter_transport_count.x_axis = _x_axis                                         
                    if re.search(r'emitter_emitter_transport_count.yml', file) and 'all_data_points' not in root:
                        sec_emitter_emitter_transport_count = sec_exciton_molpairs.m_create(Emitter_emitter_transport_count)
                        file_emitter_emitter_transport_count = yaml.safe_load(f)
                        if 'dexter_S1S1' in file_emitter_emitter_transport_count:
                            _dexter_s1s1 = file_emitter_emitter_transport_count['dexter_S1S1']
                            sec_emitter_emitter_transport_count.dexter_s1s1 = _dexter_s1s1
                        if 'dexter_S1T1' in file_emitter_emitter_transport_count:
                            _dexter_s1t1 = file_emitter_emitter_transport_count['dexter_S1T1']
                            sec_emitter_emitter_transport_count.dexter_s1t1 = _dexter_s1t1
                        if 'dexter_T1S1' in file_emitter_emitter_transport_count:
                            _dexter_t1s1 = file_emitter_emitter_transport_count['dexter_T1S1']
                            sec_emitter_emitter_transport_count.dexter_t1s1 = _dexter_t1s1
                        if 'dexter_T1T1' in file_emitter_emitter_transport_count:
                            _dexter_t1t1 = file_emitter_emitter_transport_count['dexter_T1T1']
                            sec_emitter_emitter_transport_count.dexter_t1t1 = _dexter_t1t1
                        if 'foerster_S1S1' in file_emitter_emitter_transport_count:
                            _foerster_s1s1 = file_emitter_emitter_transport_count['foerster_S1S1']
                            sec_emitter_emitter_transport_count.foerster_s1s1 = _foerster_s1s1
                        if 'foerster_S1T1' in file_emitter_emitter_transport_count:
                            _foerster_s1t1 = file_emitter_emitter_transport_count['foerster_S1T1']
                            sec_emitter_emitter_transport_count.foerster_s1t1 = _foerster_s1t1
                        if 'foerster_T1S1' in file_emitter_emitter_transport_count:
                            _foerster_t1s1 = file_emitter_emitter_transport_count['foerster_T1S1']
                            sec_emitter_emitter_transport_count.foerster_t1s1 = _foerster_t1s1
                        if 'foerster_T1T1' in file_emitter_emitter_transport_count:
                            _foerster_t1t1 = file_emitter_emitter_transport_count['foerster_T1T1']
                            sec_emitter_emitter_transport_count.foerster_t1t1 = _foerster_t1t1
                        if 'x_axis' in file_emitter_emitter_transport_count:
                            _x_axis = file_emitter_emitter_transport_count['x_axis']
                            sec_emitter_emitter_transport_count.x_axis = _x_axis
                    if re.search(r'host_host_transport_count.yml', file) and 'all_data_points' not in root:
                        sec_host_host_transport_count = sec_exciton_molpairs.m_create(Host_host_transport_count)
                        file_host_host_transport_count = yaml.safe_load(f)
                        if 'dexter_S1S1' in file_host_host_transport_count:
                            _dexter_s1s1 = file_host_host_transport_count['dexter_S1S1']
                            sec_host_host_transport_count.dexter_s1s1 = _dexter_s1s1 
                        if 'dexter_S1T1' in file_host_host_transport_count:
                            _dexter_s1t1 = file_host_host_transport_count['dexter_S1T1']
                            sec_host_host_transport_count.dexter_s1t1 = _dexter_s1t1
                        if 'dexter_T1S1' in file_host_host_transport_count:
                            _dexter_t1s1 = file_host_host_transport_count['dexter_T1S1']
                            sec_host_host_transport_count.dexter_t1s1 = _dexter_t1s1
                        if 'dexter_T1T1' in file_host_host_transport_count:
                            _dexter_t1t1 = file_host_host_transport_count['dexter_T1T1']
                            sec_host_host_transport_count.dexter_t1t1 = _dexter_t1t1
                        if 'foerster_S1S1' in file_host_host_transport_count:
                            _foerster_s1s1 = file_host_host_transport_count['foerster_S1S1']
                            sec_host_host_transport_count.foerster_s1s1 = _foerster_s1s1
                        if 'foerster_S1T1' in file_host_host_transport_count:
                            _foerster_s1t1 = file_host_host_transport_count['foerster_S1T1']
                            sec_host_host_transport_count.foerster_s1t1 = _foerster_s1t1
                        if 'foerster_T1S1' in file_host_host_transport_count:
                            _foerster_t1s1 = file_host_host_transport_count['foerster_T1S1']
                            sec_host_host_transport_count.foerster_t1s1 = _foerster_t1s1
                        if 'foerster_T1T1' in file_host_host_transport_count:
                            _foerster_t1t1 = file_host_host_transport_count['foerster_T1T1']
                            sec_host_host_transport_count.foerster_t1t1 = _foerster_t1t1
                        if 'x_axis' in file_host_host_transport_count:
                            _x_axis = file_host_host_transport_count['x_axis']
                            sec_host_host_transport_count.x_axis = _x_axis
                
                if 'runtime_analysis' in root:
                    
                    if not runtime_analysis_hasrun:    
                        sec_runtime_analysis = sec_experiments.m_create(Runtime_analysis)
                        runtime_analysis_hasrun = True
                    if re.search(r'event_counts_by_type_\d+', file) and not 'all_data_points' in root:
                        sec_event_counts_by_type = sec_runtime_analysis.m_create(Event_counts_by_type)
                        for i, line in enumerate(f):
                            line = line.lower()
                            parts = line.split(': ')
                            if 'dexter eeq' in line:
                                sec_event_counts_by_type.dexter_eeq = float(parts[1])
                            if 'dexter ept' in line:
                                sec_event_counts_by_type.dexter_ept = float(parts[1])
                            if 'spq' in line:
                                sec_event_counts_by_type.spq = float(parts[1])       
                            if 'sta' in line:
                                sec_event_counts_by_type.sta = float(parts[1])
                            if 'tpq' in line:
                                sec_event_counts_by_type.tpq = float(parts[1])
                            if 'tta' in line:
                                sec_event_counts_by_type.tta = float(parts[1])
                            if 'ttf' in line:
                                sec_event_counts_by_type.ttf = float(parts[1])
                            if 'eject chg' in line:
                                sec_event_counts_by_type.eject_chg = float(parts[1])
                            if 'inject e' in line:
                                sec_event_counts_by_type.inject_e = float(parts[1])
                            if 'inject h' in line:
                                sec_event_counts_by_type.inject_h = float(parts[1])
                            if 'move chg' in line:
                                sec_event_counts_by_type.move_chg = float(parts[1])
                            if 'move exc dexter' in line:
                                sec_event_counts_by_type.move_exc_dexter = float(parts[1])
                            if 'move exc foerster' in line:
                                sec_event_counts_by_type.move_exc_foerster = float(parts[1])
                            if 'move+flip exc dexter' in line:
                                sec_event_counts_by_type.move_flip_exc_dexter = float(parts[1])
                            if 'move+flip exc foerster' in line:
                                sec_event_counts_by_type.move_flip_exc_foerster = float(parts[1])
                            if 'prtclrst' in line:
                                sec_event_counts_by_type.prtclRst = float(parts[1])
                            if 'rad decay' in line:
                                sec_event_counts_by_type.rad_decay = float(parts[1])
                            if 'recombination s1' in line:
                                sec_event_counts_by_type.recombination_s1 = float(parts[1])
                            if 'recombination t1' in line:
                                sec_event_counts_by_type.recombination_t1 = float(parts[1])
                            if 'seperate eh' in line:
                                sec_event_counts_by_type.seperate_eh = float(parts[1])
                            if 'seperate he' in line:
                                sec_event_counts_by_type.seperate_he = float(parts[1])
                            if 'setmult' in line:
                                sec_event_counts_by_type.setMult = float(parts[1])
                            if 'spin flip exc' in line:
                                sec_event_counts_by_type.spin_flip_exc = float(parts[1])
                            if 'thermal_decay' in line:
                                sec_event_counts_by_type.thermal_decay = float(parts[1])
                if 'device_data' in root:
                    if not device_data_hasrun:
                        sec_device_data = sec_material.m_create(Device_data)
                        device_data_hasrun = True
                    
                    if re.search(r'coord_\d+', file) and 'all_data_points' not in root:
                        sec_coordinates = sec_device_data.m_create(Coordinates)
                        _coordinates = []
                        for i, line in enumerate(f):
                            parts = line.split()
                            parts = [float(p) for p in parts]
                            _coordinates.append(parts) 
                        sec_coordinates.coordinates = _coordinates      


class LightforgeParser():

    def parse(self, filepath, archive, logger):
        sec_program = archive.m_setdefault('run.program')
        sec_program.name = "Lightforge"

        sec_workflow = archive.m_create(Workflow)
        sec_workflow.type = 'single_point'
        mainfile = Path(filepath)
        
        
        
        DetailedParser(mainfile, archive)
