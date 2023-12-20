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
    text = Quantity(type=str)

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

class Foerster_expansion_errors(MSection):
    m_def = Section(validate=False)
    s1s1_0_0 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    t1t1_0_0 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    s1t1_0_0 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    t1s1_0_0 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    s1s1_0_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    t1t1_0_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    s1t1_0_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    t1s1_0_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    s1s1_1_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    t1t1_1_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    s1t1_1_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')
    t1s1_1_1 = Quantity(type=np.float64, description='cumulative error of Gaussian Foerster rate expansion')

class Dexter_and_foerster(MSection):
    m_def = Section(validate=False)
    name = Quantity(type=str, description='name of the file, which indicates the type of the "values"-matrix.')
    values = Quantity(type=np.float64, shape=['*', 2], description='1st row is rate in 1/s, 2nd row is pair distance in nm.')

class Foerster(MSection):
    m_def = Section(validate=False)
    dexter_and_foerster = SubSection(sub_section=Dexter_and_foerster.m_def, repeats=True)
    foerster_expansion_errors = SubSection(sub_section=Foerster_expansion_errors.m_def, repeats=False)
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

class Layers(MSection):
    m_def = Section(validate=False)

    lf_thickness = Quantity(type=np.float64)
    morphology_input_mode = Quantity(type=str)
    lf_material = Quantity(type=str, shape=['*', '*'], description='''Indicates what type of material this layer is made of (host, emitter etc) and its concentration 
                                                                    with 1.0 being max concentration.''')
class Pair_input(MSection):
    m_def = Section(validate=False)
    
    molecule_1_type = Quantity(type=str)
    molecule_2_type = Quantity(type=str)
    hole_transfer_integrals_wf_decay_length = Quantity(type=np.float64)
    hole_transfer_ingetrals_maximum_ti = Quantity(type=np.float64)
    electron_transfer_integrals_wf_decay_length = Quantity(type=np.float64)
    electron_transfer_integrals_maximum_ti = Quantity(type=np.float64)
    dexter_transfer_integrals_wf_decay_length = Quantity(type=np.float64)
    dexter_transfer_integrals_maximum_ti = Quantity(type=np.float64)
    
    lf_QP_output = Quantity(type=str)
class LF_materials(MSection):
    m_def = Section(validate=False)
    
    name = Quantity(type=str, description='Name of the material.')
    lf_exciton_preset = Quantity(type=str)
    lf_energies = Quantity(type=np.float64, shape=[2], description='1st value is HOMO in eV, 2nd one is LUMO in eV.')
    lf_sigma = Quantity(type=np.float64, shape=[2], description='1st value is HOMO disorder in eV, 2nd value is LUMO disorder in eV.')
    lf_lambda = Quantity(type=np.float64, shape=[2], description='1st value is HOMO reorganization energy in eV, 2nd value is LUMO reorganziation energy in eV.')
class Settings(MSection):
    m_def = Section(validate=False)
    
    pbc = Quantity(type=str)
    excitonics = Quantity(type=str)
    connect_electrodes = Quantity(type=str)
    coulomb_mesh = Quantity(type=str)
    particles_holes = Quantity(type=str)
    particles_electrons = Quantity(type=str)
    particles_excitons = Quantity(type=str)
    morphology_width = Quantity(type=np.float64)
    input_mode_transport = Quantity(type=str, description='''QP = QuantumPatch, PAR = parametric, eaip = Electron Affinity/Ionozation Potential, sig = sigma = energy disorder,
                                                            l = lambda = reorganziation energy.''')
    lf_neighbours = Quantity(type=np.float64)
    transfer_integral_source = Quantity(type=str)
    electrode_workfunction = Quantity(type=np.float64, shape=[2], description='''1st value is workfunction in eV of the electrode that is attached before the layer
                                                                                 that was defined first (defined under settings.layers.lf_material). 2nd electrode
                                                                                 is attached after the last layer.''')
    electrode_coupling_model = Quantity(type=str, shape=[2])
    electrode_wf_decay_length = Quantity(type=np.float64, shape=[2])
    electrode_coupling = Quantity(type=np.float64, shape=[2])
    lf_simulations = Quantity(type=int, description='number of simulation runs.')
    lf_measurement = Quantity(type=str)
    lf_temperature = Quantity(type=np.float64, shape=['*'])
    lf_field_direction = Quantity(type=np.float64, shape=['*'])
    lf_field_strength = Quantity(type=np.float64, shape=['*'], description='array for applied electric fields in V/nm.')
    lf_initial_holes = Quantity(type=int)
    lf_initial_electrons = Quantity(type=int)
    lf_iv_fluctuation = Quantity(type=np.float64, description='iv-fluctuation as convergence criterion.')
    lf_max_iterations = Quantity(type=int, description='max number of iterations as convergence criterion.')
    
    lf_ti_prune = Quantity(type=str)
    lf_expansion_scheme = Quantity(type=str)
    lf_qp_output_files_name = Quantity(type=str, description='files from Nanomatch GmbH software "QuantumPatch" (qp) are used.')
    lf_qp_output_files_output_zip = Quantity(type=str)
    lf_molecule_pdb = Quantity(type=str, description='name of file that provides data about the molecules.')
    qp_output_sigma = Quantity(type=str)
    qp_output_eaip = Quantity(type=str)
    lf_rates = Quantity(type=str)
    lf_superexchange = Quantity(type=str)
    lf_epsilon_material = Quantity(type=np.float64)

    pair_input = SubSection(sub_section=Pair_input.m_def, repeats=True)
    layers = SubSection(sub_section=Layers.m_def, repeats=True)
    materials = SubSection(sub_section=LF_materials.m_def, repeats=True)

class Run_lf_slr(MSection):
    m_def = Section(validate=False)

class Files_for_kmc(MSection):
    m_def = validate(Section=False)


class Input(MSection):
    m_def = Section(validate=False)
    
    settings = SubSection(sub_section=Settings.m_def, repeats=False)
    run_lf_slr = SubSection(sub_section=Run_lf_slr.m_def, repeats=False)
    files_for_kmc = SubSection(sub_section=Files_for_kmc.m_def, repeats=False)

class LF_add_info(MSection):
    m_def = Section(validate=False)

    lf_layer_id = Quantity(type=np.float64)
    n_layer_sites = Quantity(type=np.float64)
    sites_end_idx_in_device = Quantity(type=np.float64)
    sites_start_idx_in_device = Quantity(type=np.float64)
    add_info_thickness = Quantity(type=np.float64)
    add_info_x_boundaries = Quantity(type=np.float64, shape=[2])
class Material_data(MSection):
    m_def = Section(validate=False)
    
    lf_add_info = SubSection(sub_section=LF_add_info.m_def, repeats=True)
class Lightforge_data(MSection):
    m_def = Section(validate=False)
    
    material_data = SubSection(sub_section=Material_data.m_def, repeats=False)
class LightforgeCalculation(simulation.calculation.Calculation):
    m_def = Section(validate=False, extends_base_section=True)    
    
    experiments = SubSection(sub_section=Experiments.m_def, repeats=False)
    material = SubSection(sub_section=Material.m_def, repeats=False)
    input = SubSection(sub_section=Input.m_def, repeats=False)
    lightforge_data = SubSection(sub_section=Lightforge_data.m_def, repeats=False)

m_package.__init_metainfo__()
