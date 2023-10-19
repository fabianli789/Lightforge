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
from nomad.datamodel import results, optimade
#from nomad.datamodel.metainfo.workflow import Workflow


m_package = Package()
class Mobility(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=['*'], description='The number of (repeating) subsections of "mobility" corresponds to the number of applied electric fields, which are given in the "run.calculation.input" section, with subsection "0" being the first electric field.' )
    mobilities_all_fields = Quantity(type=np.float64, shape=['*', 3], description='1st column is electric field^0.5 in (V/cm)^0.5, 2nd column is mobility in cm^2/(Vs).')
class IV(MSection):
    m_def =  Section(validate=False)
    name = Quantity(type=str)
    iv_all_fields = Quantity(type=np.float64, shape=['*,3'], description='Voltage (in Volts) over current density (in mA/cm2).1st column is voltage, 3rd column is current density.')
    iv_all_fields_log  = Quantity(type=np.float64, shape=['*, 2'])
class IQE2(MSection):
    m_def = Section(validate=False)
    name = Quantity(type=str)
    iqe2_all_fields = Quantity(type=np.float64, shape=['*', 2], description='1st column is voltage in Volts, 2nd column is IQE.')
    iqe2_all_currents = Quantity(type=np.float64, shape=['*',2], description='1st column is current density in mA/cm2, 2nd column is the Internal Quantum Efficiency (IQE). The current densities are a factor 10 too large.')
class Current_density(MSection):
    m_def = Section(validate=False)
    value = Quantity(type=np.float64, shape=['*'],  description='The number of (repeating) subsections of "current_density" corresponds to the number of applied electric fields, which are given in the "run.calculation.input" section, with subsection "0" being the first electric field.')

class Current_characteristics(MSection):
    m_def = Section(validate=False)
    current_density = SubSection(sub_section=Current_density.m_def, repeats=True)    
    IQE2 = SubSection(sub_section=IQE2.m_def, repeats=False)
    IV = SubSection(sub_section=IV.m_def, repeats=False)
    mobility = SubSection(sub_section=Mobility.m_def, repeats=True)
class Experiments(MSection):
    m_def = Section(validate=False)
    current_characteristics = SubSection(sub_section=Current_characteristics.m_def, repeats=False)
class Material(MSection):
    m_def = Section(validate=False)
    
class Input(MSection):
    m_def = Section(validate=False)

class Calculation(simulation.calculation.Calculation):
    m_def = Section(validate=False, extends_base_section=True)    
    experiments = SubSection(sub_section=Experiments.m_def, repeats=False)
    material = SubSection(sub_section=Material.m_def, repeats=False)
    input = SubSection(sub_section=Input.m_def, repeats=False)
    
class Run(simulation.run.Run):
    m_def = Section(validate=False, extends_base_section=True)

m_package.__init_metainfo__()
