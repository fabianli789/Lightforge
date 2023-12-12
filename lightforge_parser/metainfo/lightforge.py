#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np
from nomad.metainfo import MSection, Section, SubSection, Quantity, Package
from nomad.datamodel.metainfo import simulation
from nomad.datamodel import results
#from nomad.datamodel.metainfo.workflow import Workflow


m_package = Package()

class Site_energies(MSection):
    m_def = Section(validate=False)
    site_energies = Quantity(type=np.float64, shape=['*', 2])

class Mol_types(MSection):
    m_def = Section(validate=False)
    mol_types = Quantity(type=np.float64, shape=['*'])

class Coordinates(MSection):
    m_def = Section(validate=False)
    coordinates = Quantity(type=np.float64, shape=['*', 3])


class Device_data(MSection):
    m_def = Section(validate=False)
    coordinates = SubSection(sub_section=Coordinates.m_def, repeats=True)
    mol_types = SubSection(sub_section=Mol_types.m_def, repeats=True)
    site_energies = SubSection(sub_section=Site_energies.m_def, repeats=True)
class Electrodes(MSection):
    m_def = Section(validate=False)

class Energy_levels(MSection):
    m_def = Section(validate=False)

class Exciton_separation(MSection):
    m_def = Section(validate=False)

class Foerster(MSection):
    m_def = Section(validate=False)
    
class Emitter_emitter_transport_count(MSection):
    m_def = Section(validate=False)
    dexter_s1s1 = Quantity(type=np.float64, shape=['*'])
    dexter_s1t1 = Quantity(type=np.float64, shape=['*'])
    dexter_t1s1 = Quantity(type=np.float64, shape=['*'])
    dexter_t1t1 = Quantity(type=np.float64, shape=['*'])
    foerster_s1s1 = Quantity(type=np.float64, shape=['*'])
    foerster_s1t1 = Quantity(type=np.float64, shape=['*'])
    foerster_t1s1 = Quantity(type=np.float64, shape=['*'])
    foerster_t1t1 = Quantity(type=np.float64, shape=['*'])
    x_axis = Quantity(type=np.float64, shape=['*'])
class Host_host_transport_count(Emitter_emitter_transport_count):
    m_def = Section(validate=False)

class Host_emitter_transport_count(Emitter_emitter_transport_count):
    m_def = Section(validate=False)

class Exciton_molpairs(MSection):
    m_def = Section(validate=False)
    emitter_emitter_transport_count = SubSection(sub_section=Emitter_emitter_transport_count.m_def, repeats=False)
    host_emitter_transport_count = SubSection(sub_section=Host_emitter_transport_count.m_def, repeats=False)  # might be wrong
    host_host_transport_count = SubSection(sub_section=Host_host_transport_count.m_def, repeats=False)
class Quenching_density_average(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=[4, '*'], description="""1st row is device length in nm,
                                                                     2nd row is average number of quenched
                                                                     excitons, 3rd row is average 
                                                                     number of electrons, 4th row is 
                                                                     average number of holes.""")

class Photon_creation_density_average(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=[3, '*'], description="""1st row is device length in nm, 
                                                                    2nd row is average number of emitted photons, 
                                                                    3rd row is average number of excitons created. 
                                                                    1st repeating subsection corresponds to the 1st 
                                                                    electric field applied etc.""")

class Exciton_decay_density_average(MSection):
    m_def = Section(validate=False)
    
    cs = Quantity(type=np.float64, shape=['*'], description='avg. number of charge separations.')
    eet = Quantity(type=np.float64, shape=['*'])
    ept = Quantity(type=np.float64, shape=['*'])
    ptq = Quantity(type=np.float64, shape=['*'], description='polaron-triplet quenching.')
    spq = Quantity(type=np.float64, shape=['*'], description='singlet-polaron quenching.')
    ssa = Quantity(type=np.float64, shape=['*'], description='singlet-singlet annihilation')        
    sta = Quantity(type=np.float64, shape=['*'], description='singlet-triplet annihilation')
    tpq = Quantity(type=np.float64, shape=['*'], description='triplet-polaron quenching')
    tsa = Quantity(type=np.float64, shape=['*'], description='triplet-singlet annihilation')
    tta = Quantity(type=np.float64, shape=['*'], description='triplet-triplet annihilation')
    ttf = Quantity(type=np.float64, shape=['*'], description='triplet-triplet fusion')
    radiative = Quantity(type=np.float64, shape=['*'], description='avg. no. of radiatively decayed photons')
    thermal = Quantity(type=np.float64, shape=['*'], description='thermal quenching')    
    photon = Quantity(type=np.float64, shape=['*'], description='number of photons created')
    recombination = Quantity(type=np.float64, shape=['*'], description='avg. number of recombinations')
    x_axis = Quantity(type=np.float64, shape=['*'], description='device length in nanometers.')
    
class Charge_density_average(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=[3, '*'], description="""1st row is device length in nm, 
                                                                     2nd row is average number of electrons, 
                                                                     3rd row is average number of holes. 
                                                                     1st repeating subsection corresponds to 
                                                                     the first electric field applied etc.""")

class Mobility(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=['*'], description="""The number of (repeating) subsections of 
                                                                  "mobility" corresponds to the number of 
                                                                  applied electric fields, which are given 
                                                                  in the "run.calculation.input" section, 
                                                                  with subsection "0" being the first 
                                                                  electric field.""")
    mobilities_all_fields = Quantity(type=np.float64, shape=['*', 3], description="""1st column is 
                                                                                     electric field^0.5 
                                                                                     in (V/cm)^0.5, 2nd 
                                                                                     column is mobility 
                                                                                     in cm^2/(Vs).""")
class IV(MSection):
    m_def =  Section(validate=False)
    name = Quantity(type=str)
    iv_all_fields = Quantity(type=np.float64, shape=['*', 3], description="""Voltage (in Volts) over 
                                                                            current density (in mA/cm2).
                                                                            1st column is voltage, 3rd 
                                                                            column is current density.""")
    iv_all_fields_log  = Quantity(type=np.float64, shape=['*, 2'])
class IQE2(MSection):
    m_def = Section(validate=False)
    name = Quantity(type=str)
    iqe2_all_fields = Quantity(type=np.float64, shape=['*', 2], description="""1st column is voltage in 
                                                                               Volts, 2nd column is IQE.""")
    iqe2_all_currents = Quantity(type=np.float64, shape=['*',2], description="""1st column is current 
                                                                                density in mA/cm2, 2nd 
                                                                                column is the Internal 
                                                                                Quantum Efficiency (IQE).
                                                                                The current densities are 
                                                                                a factor 10 too large.""")
class Current_density(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=['*'],  description="""The number of (repeating) 
                                                                   subsections of "current_density" 
                                                                   corresponds to the number of applied 
                                                                   electric fields, which are given in 
                                                                   the "run.calculation.input" section, 
                                                                   with subsection "0" being the first 
                                                                   electric field.""")

class Event_counts_by_type(Exciton_decay_density_average):
    m_def = Section(validate=False)
    dexter_eeq = Quantity(type=np.float64)
    dexter_ept = Quantity(type=np.float64)
    eject_chg = Quantity(type=np.float64)
    inject_e = Quantity(type=np.float64)
    inject_h = Quantity(type=np.float64)
    move_chg = Quantity(type=np.float64)
    move_exc_dexter = Quantity(type=np.float64)
    move_exc_foerster = Quantity(type=np.float64)
    move_flip_exc_dexter = Quantity(type=np.float64)
    move_flip_exc_foerster = Quantity(type=np.float64)
    prtclRst = Quantity(type=np.float64)
    rad_decay = Quantity(type=np.float64)
    recombination_s1 = Quantity(type=np.float64)
    recombination_t1 = Quantity(type=np.float64)
    seperate_eh = Quantity(type=np.float64)
    seperate_he = Quantity(type=np.float64)
    setMult = Quantity(type=np.float64)
    spin_flip_exc = Quantity(type=np.float64)
    thermal_decay = Quantity(type=np.float64)
class Particle_densities(MSection):
    m_def = Section(validate=False)
    charge_density_average = SubSection(sub_section=Charge_density_average.m_def, repeats=True)
    exciton_decay_density_average = SubSection(sub_section=Exciton_decay_density_average.m_def, repeats=True)
    photon_creation_density_average = SubSection(sub_section=Photon_creation_density_average.m_def, repeats=True)
    quenching_density_average = SubSection(sub_section=Quenching_density_average.m_def, repeats=True)
    exciton_molpairs = SubSection(sub_section=Exciton_molpairs.m_def, repeats=False)
class Current_characteristics(MSection):
    m_def = Section(validate=False)
    current_density = SubSection(sub_section=Current_density.m_def, repeats=True)    
    IQE2 = SubSection(sub_section=IQE2.m_def, repeats=False)
    IV = SubSection(sub_section=IV.m_def, repeats=False)
    mobility = SubSection(sub_section=Mobility.m_def, repeats=True)

class Runtime_analysis(MSection):
    m_def = Section(validate=False)
    event_counts_by_type = SubSection(sub_section=Event_counts_by_type.m_def, repeats=True)

class Experiments(MSection):
    m_def = Section(validate=False)
    current_characteristics = SubSection(sub_section=Current_characteristics.m_def, repeats=False)
    particle_densities = SubSection(sub_section=Particle_densities.m_def, repeats=False)
    runtime_analysis = SubSection(sub_section=Runtime_analysis.m_def, repeats=False)
class Material(MSection):
    m_def = Section(validate=False)
    device_data = SubSection(sub_section=Device_data.m_def, repeats=False)
    electrodes = SubSection(sub_section=Electrodes.m_def, repeats=False)
    energy_levels = SubSection(sub_section=Energy_levels.m_def, repeats=False)
    exciton_separation = SubSection(sub_section=Exciton_separation.m_def, repeats=False)
    foerster = SubSection(sub_section=Foerster.m_def, repeats=False)
class Input(MSection):
    m_def = Section(validate=False)

class LightforgeCalculation(simulation.calculation.Calculation):
    m_def = Section(validate=False, extends_base_section=True)    
    experiments = SubSection(sub_section=Experiments.m_def, repeats=False)
    material = SubSection(sub_section=Material.m_def, repeats=False)
    input = SubSection(sub_section=Input.m_def, repeats=False)


m_package.__init_metainfo__()
