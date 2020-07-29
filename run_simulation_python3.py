#!/usr/bin/python

from neuron import h
from neuron import gui
import numpy as np
import os
import pickle
import pickle
from scipy.stats.stats import pearsonr
import sys
sys.path.append(os.getcwd())
import itertools

print ('Setting cell model...')

neuron_num = 9
synaptic_map_file_name = 'synapse_location_for_neuron_{}.pickle'.format(neuron_num)
rate_E = 1.75
rate_I = 10

h.load_file('model_files/{}.hoc'.format(neuron_num))
h.v_init = -80
h.dt = 0.025

def set_cell(neuron_num):
  h.load_file('model_files/{}.hoc'.format(neuron_num))
  h('''
  objref somatic, apical, axonal, basal
  somatic = new SectionList()
  basal = new SectionList()
  apical = new SectionList()
  axonal = new SectionList()
  celsius=37

  proc distribute(){local x localobj sl
          strdef stmp,distfunc,mech
          sl = $o1
          mech = $s2
          distfunc = $s3
          sprint(distfunc,"%%s %s(%%f) = %s",mech, distfunc)
          forsec sl for(x) {
              sprint (stmp,distfunc,secname(),x,distance(x))
              execute(stmp)
          }
      }

  create axon[2]''')
  
  
  h('''
    access axon[0] 
    {
      L = 30                                                              
      diam = 0.66
      nseg = 30                                         
      all.append()                                                            
      axonal.append()                                                         
    }
    
    access axon[1] 
    {                                                            
        L = 30                                                                  
        diam = 0.33
        nseg = 30                                     
        all.append()                                                            
        axonal.append()                                                         
    }                                                                           
    nSecAxonal = 2                                                              
        
    soma connect axon[0](0), 1                                           
    axon[0] connect axon[1](0), 1
  ''')

  h("celsius = 37.0")
  for sec in h.all:
    if sec.name().startswith('soma') or sec.name().startswith('axon'): continue
    h("{}.nseg = {}".format(sec.name(), max(eval('h.{}.L'.format(sec.name())), 1)))

  for sec in h.all:
    if sec.name().startswith('dend_5') or sec.name().startswith('apic'):
      h('access {}'.format(sec.name()))
      h('apical.append()')
      h.pop_section()
      
    elif sec.name().startswith('dend_7'):
      h('access {}'.format(sec.name()))
      h('basal.append()')
      h.pop_section()
      
    elif 'soma' in sec.name() or sec.name().startswith('dend_0'):
      h('access {}'.format(sec.name()))
      h('somatic.append()')
      h.pop_section()
      
    elif 'axon' in sec.name():
      h('access {}'.format(sec.name()))
      h('axonal.append()')
      h.pop_section()
      
    else:
      print('no section list for {}'.format(sec.name()))
    
  h('''
    forsec all {
      insert pas
    }

    forsec all {
      e_pas = -75 
    }

    forsec all {
      Ra = 150
    }

    forsec all {
      cm = 1
    }

    forsec all {
      g_pas = 1.0 / 12000.0
    }

    forsec apical {
      cm = 1
    }

    forsec basal {
      cm = 1
    }
  ''')

def active_soma():
  h('''
    forsec somatic {
      insert Ca_HVA
      insert SKv3_1
      insert SK_E2
      insert Ca_LVAst
      insert Ih
      insert NaTs2_t
      insert CaDynamics_E2
    }

    forsec axonal {
      insert Ca_HVA
      insert SKv3_1
      insert SK_E2
      insert CaDynamics_E2
      insert Nap_Et2
      insert K_Pst
      insert K_Tst
      insert Ca_LVAst
      insert NaTa_t
      
    }

    forsec axonal {
      ena = 50 
      ek = -85 
    }

    forsec somatic {
      ena = 50 
    }

    forsec somatic {
      ek = -85 
    }
    
    access soma
    distance()

    distribute(axonal,"gNaTa_tbar_NaTa_t","(0.0 * %g + 1.0)*3.429725")
    distribute(axonal,"gK_Tstbar_K_Tst","(0.0 * %g + 1.0)*0.001035")
    distribute(axonal,"gamma_CaDynamics_E2","(0.0 * %g + 1.0)*0.016713")
    distribute(axonal,"gNap_Et2bar_Nap_Et2","(0.0 * %g + 1.0)*0.009803")
    distribute(axonal,"gSK_E2bar_SK_E2","(0.0 * %g + 1.0)*0.008085")
    distribute(axonal,"gCa_HVAbar_Ca_HVA","(0.0 * %g + 1.0)*0.000306")
    distribute(axonal,"gK_Pstbar_K_Pst","(0.0 * %g + 1.0)*0.959296")
    distribute(axonal,"gSKv3_1bar_SKv3_1","(0.0 * %g + 1.0)*0.094971")
    distribute(axonal,"decay_CaDynamics_E2","(0.0 * %g + 1.0)*384.114655")
    distribute(axonal,"gCa_LVAstbar_Ca_LVAst","(0.0 * %g + 1.0)*0.000050")
    distribute(somatic,"gamma_CaDynamics_E2","(0.0 * %g + 1.0)*0.000533")
    distribute(somatic,"gSKv3_1bar_SKv3_1","(0.0 * %g + 1.0)*0.102517")
    distribute(somatic,"gSK_E2bar_SK_E2","(0.0 * %g + 1.0)*0.099433")
    distribute(somatic,"gCa_HVAbar_Ca_HVA","(0.0 * %g + 1.0)*0.000374")
    distribute(somatic,"gNaTs2_tbar_NaTs2_t","(0.0 * %g + 1.0)*0.926705")
    distribute(somatic,"gIhbar_Ih","(0.0 * %g + 1.0)*0.000080")
    distribute(somatic,"decay_CaDynamics_E2","(0.0 * %g + 1.0)*342.544232")
    distribute(somatic,"gCa_LVAstbar_Ca_LVAst","(0.0 * %g + 1.0)*0.000778")    
  ''')

def active_dendrites_tuned():
  h('''
    forsec basal {
      insert Ih
    }
    
    forsec apical {
      insert Im
      insert NaTs2_t
      insert SKv3_1
      insert Ih
    }
    
    forsec apical {
      ena = 50 
      ek = -85 
    }
    
    access soma
    distance()

    distribute(basal,"gIhbar_Ih","(0.0 * %g + 1.0)*0.000080")
    distribute(apical,"gNaTs2_tbar_NaTs2_t","(0.0 * %g + 1.0)*0.008009") // based on ttx comparison
    // distribute(apical,"gNaTs2_tbar_NaTs2_t","(0.0 * %g + 1.0)*0.012009") // original
    distribute(apical,"gSKv3_1bar_SKv3_1","(0.0 * %g + 1.0)*0.000513")
    distribute(apical,"gIhbar_Ih","(-0.869600 + 2.087000*exp((%g-0.000000)*0.020100))*0.000080") // based on attenuation comparison
    // distribute(apical,"gIhbar_Ih","(-0.869600 + 2.087000*exp((%g-0.000000)*0.003100))*0.000080") // original
    distribute(apical,"gImbar_Im","(0.0 * %g + 1.0)*0.000740")
  ''')

set_cell(neuron_num)
active_soma()
active_dendrites_tuned()

print ('Setting parameters...')

duration_of_simulation = 1000

number_of_inhibitory_synapses_on_soma = 100
(spine_locations, full_IS_locations) = pickle.load(open('model_files/{}'.format(synaptic_map_file_name), 'rb'), encoding='latin1')

full_IS_locations = list(full_IS_locations)

for loc in np.linspace(0, 1, number_of_inhibitory_synapses_on_soma):
  full_IS_locations.append(('soma', loc, 'soma',0,0))      

h.tstop = duration_of_simulation
number_of_E_synapses = len(spine_locations)
number_of_I_synapses = len(full_IS_locations)

print ('Set recording...')

voltage_vectors = []

print ('Placing synapses...')

background_E_events = np.random.poisson(rate_E / 1000.0, size=(number_of_E_synapses, duration_of_simulation))
background_I_events = np.random.poisson(rate_I / 1000.0, size=(number_of_I_synapses, duration_of_simulation))

eSynlist = []
eNetconlist = []
iSynlist = []
iNetconlist = []
E_vcs = []
I_vcs = []
E_vcs_events = []
I_vcs_events = []

def placeNMDA(location, conductance):
  eSynlist.append(h.ProbAMPANMDA2_RATIO(float(location)))
    
  eSynlist[-1].gmax = conductance
  eSynlist[-1].mgVoltageCoeff = 0.08
  eNetconlist.append(h.NetCon(E_vcs[-1], eSynlist[-1]))
  eNetconlist[-1].weight[0] = 1
  eNetconlist[-1].delay = 0

def placeGABA(location, conductance):
  iSynlist.append(h.ProbUDFsyn2_lark(float(location)))
  iSynlist[-1].tau_r = 0.18
  iSynlist[-1].tau_d = 5
  iSynlist[-1].e = - 80
  iSynlist[-1].gmax = conductance
  iNetconlist.append(h.NetCon(I_vcs[-1], iSynlist[-1]))
  iNetconlist[-1].weight[0] = 1
  iNetconlist[-1].delay = 0

for sec in h.all:
  if sec.name().startswith('soma') or sec.name().startswith('axon'): continue
  h("{}.nseg = {}".format(sec.name(), max(eval('h.{}.L'.format(sec.name())), 1)))

for sec in h.all:
  voltage_vectors.append(h.Vector())
  voltage_vectors[-1].record(eval('h.{}({})._ref_v'.format(sec.name(), 0.5)))

spine_ind = 0
spine_firing_times = []
for sec, loc, branch_type, spine_head_width, spine_volume, spine_length in spine_locations:
  h("access {}".format(sec))
  h("nseg = {}".format(max(h.L, 1)))
  E_vcs.append(h.VecStim())
  E_vcs_events.append(h.Vector())
  spine_firing_times.append([])
  events = np.where(background_E_events[spine_ind, :].flatten())[0] + 100
  for event in events:
    E_vcs_events[-1].append(event)
    spine_firing_times[-1].append(event)
    
  E_vcs[-1].play(E_vcs_events[-1])
  placeNMDA(loc, 0.0004)  
  spine_ind += 1
  h.pop_section()        

IS_ind = 0
IS_firing_times = []
for sec, loc, branch_type, IS_volume, on_spine in full_IS_locations:
  h("access {}".format(sec))
  h("nseg = {}".format(max(h.L, 1)))
  I_vcs.append(h.VecStim())
  I_vcs_events.append(h.Vector())
  IS_firing_times.append([])

  events = np.where(background_I_events[IS_ind, :].flatten())[0] + 100
  for event in events:
    I_vcs_events[-1].append(event)   
    IS_firing_times[-1].append(event)   

  I_vcs[-1].play(I_vcs_events[-1])
  h("access {}".format(sec))
  placeGABA(loc, 0.001)
  h.pop_section()        
  IS_ind += 1

firing_times = np.zeros(shape=(int(number_of_E_synapses), int(h.tstop) + 100))
for ind in range(len(spine_firing_times)):
  firing_times[ind, spine_firing_times[ind]] = 1

print ('Running...')

soma_voltageVector = h.Vector()
soma_voltageVector.record(h.soma(0.5)._ref_v)

timeVector = h.Vector()
timeVector.record(h._ref_t)

E_densities = []
I_densities = []

h.stdinit()
h.run()

voltage_vectors = np.array(voltage_vectors)
timeVector = np.array(timeVector)

