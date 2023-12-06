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
                            Quenching_density_average)


def DetailedParser(filepath, archive):
    sec_run = archive.m_create(Run)
    sec_calc = sec_run.m_create(Calculation)
    sec_experiments =  sec_calc.m_create(Experiments)
    sec_current_characteristics = sec_experiments.m_create(Current_characteristics)
    sec_particle_densities = sec_experiments.m_create(Particle_densities)
    
    sec_IQE2 = sec_current_characteristics.m_create(IQE2)
    sec_IV = sec_current_characteristics.m_create(IV)
    for root, dirs, files in sorted(os.walk(filepath.parent)):
        files = sorted(files)
        i = 0
        while i < len(files):
            if '.png' in files[i] or '.npz' in files[i]:
                files.remove(files[i])
            else:
                i += 1    
        
        for file in files:    
            with open(root +'/'+ file, 'rb') as f:
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

                        

class LightforgeParser():

    def parse(self, filepath, archive, logger):
        sec_program = archive.m_setdefault('run.program')
        sec_program.name = "Lightforge"

        sec_workflow = archive.m_create(Workflow)
        sec_workflow.type = 'single_point'
        mainfile = Path(filepath)
        
        
        
        DetailedParser(mainfile, archive)
